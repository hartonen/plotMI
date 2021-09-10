#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import seaborn as sns

from time import time
import argparse
import csv
import numpy as np
import gzip
import sys

import multiprocessing as mp

from helpers import getP_j, getMI_mn, getBC_mn, getHE_mn, getJS_mn

def plotMI():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #PARAMETERS
    parser.add_argument("--outdir",help="Full path to the output directory.",type=str)
    parser.add_argument("--seqs",help="Full path to the plain text input sequence file. Each sequence must be of same length and on its separate line.",type=str)
    parser.add_argument("--distance",help="Distance used to compare positional k-mer distributions. MI=mutual information (default), JS=Jensen-Shannon divergence, JS_inv=inverted Jensen-Shannon distance, BC=Bhattacharyya distance, BC_inv=inverted Bhattacharyya distance, HE=Hellinger distance.",type=str,choices=['MI','JS','JS_inv','BC','BC_inv','HE'],default='MI')
    parser.add_argument("--nproc",help="Number of parallel processes used when computing MI (default=1).",type=int,default=1)
    parser.add_argument("--figtype",help="png or pdf (default=png).",type=str,choices=['pdf','png'],default='png')
    parser.add_argument("--k",help="length of k-mer distributions used to calculate MI (default=3).",type=int,default=3)
    parser.add_argument("--v",help="Verbosity level, 0=none, 1=print info on screen (default=1).",type=int,choices=[0,1],default=1)
    parser.add_argument("--p",help="Multiplier for pseudocount mass added to k-mer count. Total pseudocount mass added is p*(number of sequences) (default=5).",type=float,default=5)
    parser.add_argument("--alphabet",help="A string containing each individual letter in the alphabet used (default=ACGT). NOTE! This is case-sensitive.",type=str,default="ACGT")
    parser.add_argument("--colorscale",help="If set to log, colormap is scaled logarithmically (default=lin, meaning linear scaling).",type=str,choices=['lin','log'],default='lin')
    parser.add_argument("--minmi",help="Set minimum value for colormap, helpful if you want to be sure that the minimum value is 0 (default=minimum value in MI matrix).",default=None,type=float)
    parser.add_argument("--step",help="Step size for axis ticks in MI-plot (default=20).",type=int,default=20)
    parser.add_argument("--save_distributions",help="If yes, save the positional and pairwise k-mer distributions and the MI contributions from each k-mer and position pair into a separate file (default=no). Note that this is a large file, possibly many GBs.",type=str,choices=['yes','no'],default='no')
    parser.add_argument("--randomized_pairs",help="EXPERIMENTAL FEATURE: File containing list of position pairs whose k-mer distributions will be randomized to flat uniform distribution according to given alphabet before computing MI. Position indices should start from 0 and be saved in tab-separated format with one pair on a single row.",type=str,default=None)
    
    args = parser.parse_args()

    #save the command used to evoke plotMI to a log file
    with open(args.outdir+"plotMI_log.txt",'wt') as logfile:
        logfile.write(" ".join(sys.argv))
        
    if args.randomized_pairs is None: randomized_pairs = []
    else:
        randomized_pairs = []
        with open(args.randomized_pairs,'rt') as infile:
            r = csv.reader(infile,delimiter='\t')
            for row in r: randomized_pairs.append((int(row[0]),int(row[1])))
    
    #read in the sequences and store them as strings
    start = time()
    seqs = []
    with open(args.seqs,'rt') as infile:
        r = csv.reader(infile)
        for row in r:
            if len(row)<1: continue
            seqs.append(row[0].upper())
    
    I = len(seqs) #number of sequences
    J = len(seqs[0]) #length of sequences
    p = args.p*I #pseudocount mass added to k-mers
    alphabet = ''.join(set(args.alphabet)) #alphabet used
    if args.v>0: print("Alphabet: "+alphabet)
    
    end = time()
    if args.v>0: print("Read in "+str(I)+" sequences of length "+str(J)+" in "+str(end-start)+" seconds.")
    
    start = time()
    #calculate single site frequencies of k-mers in parallel
    pool = mp.Pool(args.nproc)
    res = [pool.apply_async(getP_j,args=(seqs,j,I,J,args.k,p,alphabet)) for j in range(0,J-args.k+1)]
    P = [r.get() for r in res] #list containing the single site k-mer frequencies. Each distribution is saved as a dictionary where k-mer is key and its frequency is value
    
    pool.close()
    pool.join()
    pool.terminate()
    
    end = time()
    if args.v>0: print("Calculated the positional "+str(args.k)+"-mer frequencies in "+str(end-start)+" seconds.")
    
    start = time()
    if args.distance=='MI':
        #calculate MI in parallel
        if args.nproc>1: pool = mp.Pool(args.nproc)
        res = []
    
        M = 0 #number of pairwise k-mer distributions
        for m in range(0,J-2*args.k+1):
            M += 1
            for n in range(m+args.k,J-args.k+1):
                if args.nproc>1: res.append(pool.apply_async(getMI_mn,args=(seqs,m,n,I,J,P,args.k,p,alphabet,randomized_pairs)))
                else: res.append(getMI_mn(seqs,m,n,I,J,P,args.k,p,alphabet,randomized_pairs))
        if args.nproc>1:
            res = [r.get() for r in res]    

            pool.close()
            pool.join()
            pool.terminate()
    
        #build the MI matrix
        #save all individual position and k-mer contributions of MI into a file
        if args.save_distributions=='yes':
            with gzip.open(args.outdir+"MI_contributions.txt.gz",'wt') as outfile:
                w = csv.writer(outfile,delimiter='\t')
                w.writerow(['#i','j','a','b','MI_ij(a,b)','P_ij(a,b)','P_i(a)','P_j(b)'])
                MI = -1*np.ones(shape=(M,M))
                for j in range(0,len(res)):
            
                    inds = res[j][0]
                    mi = res[j][1]
                    MI[inds[1]-args.k,inds[0]] = mi
                    MI[inds[0],inds[1]-args.k] = mi
                    #saving all individual MI contributions to file
                    for kmer in res[j][2]: w.writerow([inds[0],inds[1],kmer[:args.k],kmer[args.k:],res[j][2][kmer],res[j][3][kmer],P[inds[0]][kmer[:args.k]],P[inds[1]][kmer[args.k:]]])
        else:
            MI = -1*np.ones(shape=(M,M))
            for j in range(0,len(res)):

                inds = res[j][0]
                mi = res[j][1]
                MI[inds[1]-args.k,inds[0]] = mi
                MI[inds[0],inds[1]-args.k] = mi
    elif args.distance=='JS' or args.distance=='JS_inv':
        #calculate Jensen-Shannon divergence in parallel
        pool = mp.Pool(args.nproc)
        res = []

        M = 0 #number of position pairs
        for m in range(0,J-2*args.k+1):
            M += 1
            for n in range(m+args.k,J-args.k+1):
                if args.distance=='JS': res.append(pool.apply_async(getJS_mn,args=(seqs,m,n,I,J,P,args.k,p,alphabet)))
                else: res.append(pool.apply_async(getJS_mn,args=(seqs,m,n,I,J,P,args.k,p,alphabet,True)))
        res = [r.get() for r in res]

        pool.close()
        pool.join()
        pool.terminate()
        #build the JS matrix
        if args.save_distributions=='yes':
            with gzip.open(args.outdir+"JS_contributions.txt.gz",'wt') as outfile:
                w = csv.writer(outfile,delimiter='\t')
                w.writerow(['#i','j','a','b','JS_ij(a,b)','P_i(a)','P_j(a)'])
                MI = -1*np.ones(shape=(M,M))
                for j in range(0,len(res)):

                    inds = res[j][0]
                    mi = res[j][1]
                    MI[inds[1]-args.k,inds[0]] = mi
                    MI[inds[0],inds[1]-args.k] = mi
                    #saving all individual JS contributions to file
                    for kmer in res[j][2]: w.writerow([inds[0],inds[1],kmer[:args.k],kmer[args.k:],res[j][2][kmer],P[inds[0]][kmer[:args.k]],P[inds[1]][kmer[args.k:]]\
])
        else:
            MI = -1*np.ones(shape=(M,M))
            for j in range(0,len(res)):

                inds = res[j][0]
                mi = res[j][1]
                MI[inds[1]-args.k,inds[0]] = mi
                MI[inds[0],inds[1]-args.k] = mi
    elif args.distance=='BC' or args.distance=='BC_inv':
        #calculate BC distance in parallel
        pool = mp.Pool(args.nproc)
        res = []

        M = 0 #number of position pairs
        for m in range(0,J-2*args.k+1):
            M += 1
            for n in range(m+args.k,J-args.k+1):
                if args.distance=='BC': res.append(pool.apply_async(getBC_mn,args=(seqs,m,n,I,J,P,args.k,p,alphabet)))
                else: res.append(pool.apply_async(getBC_mn,args=(seqs,m,n,I,J,P,args.k,p,alphabet,True)))
        res = [r.get() for r in res]

        pool.close()
        pool.join()
        pool.terminate()
        #build the BC matrix
        if args.save_distributions=='yes':
            with gzip.open(args.outdir+"BC_contributions.txt.gz",'wt') as outfile:
                w = csv.writer(outfile,delimiter='\t')
                w.writerow(['#i','j','a','a','BC_ij(a,a)','P_i(a)','P_j(a)'])
                MI = -1*np.ones(shape=(M,M))
                for j in range(0,len(res)):

                    inds = res[j][0]
                    mi = res[j][1]
                    MI[inds[1]-args.k,inds[0]] = mi
                    MI[inds[0],inds[1]-args.k] = mi
                    #saving all individual MI contributions to file
                    for kmer in res[j][2]: w.writerow([inds[0],inds[1],kmer[:args.k],kmer[args.k:],res[j][2][kmer],P[inds[0]][kmer[:args.k]],P[inds[1]][kmer[args.k:]]])
        else:
            MI = -1*np.ones(shape=(M,M))
            for j in range(0,len(res)):

                inds = res[j][0]
                mi = res[j][1]
                MI[inds[1]-args.k,inds[0]] = mi
                MI[inds[0],inds[1]-args.k] = mi
    elif args.distance=='HE':
        #calculate HE distance in parallel
        pool = mp.Pool(args.nproc)
        res = []

        M = 0 #number of position pairs
        for m in range(0,J-2*args.k+1):
            M += 1
            for n in range(m+args.k,J-args.k+1): res.append(pool.apply_async(getHE_mn,args=(seqs,m,n,I,J,P,args.k,p,alphabet)))
        res = [r.get() for r in res]

        pool.close()
        pool.join()
        pool.terminate()
        #build the HE matrix
        if args.save_distributions=='yes':
            with gzip.open(args.outdir+"HE_contributions.txt.gz",'wt') as outfile:
                w = csv.writer(outfile,delimiter='\t')
                w.writerow(['#i','j','a','a','HE_ij(a,a)','P_i(a)','P_j(a)'])
                MI = -1*np.ones(shape=(M,M))
                for j in range(0,len(res)):

                    inds = res[j][0]
                    mi = res[j][1]
                    MI[inds[1]-args.k,inds[0]] = mi
                    MI[inds[0],inds[1]-args.k] = mi
                    #saving all individual MI contributions to file
                    for kmer in res[j][2]: w.writerow([inds[0],inds[1],kmer[:args.k],kmer[args.k:],res[j][2][kmer],P[inds[0]][kmer[:args.k]],P[inds[1]][kmer[args.k:]]])
        else:
            MI = -1*np.ones(shape=(M,M))
            for j in range(0,len(res)):

                inds = res[j][0]
                mi = res[j][1]
                MI[inds[1]-args.k,inds[0]] = mi
                MI[inds[0],inds[1]-args.k] = mi
    end = time()
    if args.v>0: print("Computed distances between positional k-mer distributions in "+str(end-start)+" seconds.")
    
    start = time()
    #save the matrix to file
    if args.distance=='MI':
        np.savetxt(args.outdir+"MI.txt.gz",MI,delimiter='\t')
        cbar_title = 'MI'
    elif args.distance=='BC':
        np.savetxt(args.outdir+"BC.txt.gz",MI,delimiter='\t')
        cbar_title = 'Bhattacharyya distance'
    elif args.distance=='BC_inv':
        np.savetxt(args.outdir+"BC_inv.txt.gz",MI,delimiter='\t')
        cbar_title = 'inverted Bhattacharyya distance'
    elif args.distance=='HE':
        np.savetxt(args.outdir+"HE.txt.gz",MI,delimiter='\t')
        cbar_title = 'Hellinger distance'
    elif args.distance=='JS':
        np.savetxt(args.outdir+"JS.txt.gz",MI,delimiter='\t')
        cbar_title = 'Jensen-Shannon divergence'
    elif args.distance=='JS_inv':
        np.savetxt(args.outdir+"JS_inv.txt.gz",MI,delimiter='\t')
        cbar_title = 'inverted Jensen-Shannon divergence'
    
    if args.minmi!=None:
        if args.colorscale=='log': sns_plot = sns.heatmap(MI,cmap='viridis',cbar=True,cbar_kws={'label': cbar_title},vmin=args.minmi,norm=LogNorm())
        else: sns_plot = sns.heatmap(MI,cmap='viridis',cbar=True,cbar_kws={'label': cbar_title},vmin=args.minmi)
    else:
        if args.colorscale=='log': sns_plot = sns.heatmap(MI,cmap='viridis',cbar=True,cbar_kws={'label': cbar_title},norm=LogNorm())
        else: sns_plot = sns.heatmap(MI,cmap='viridis',cbar=True,cbar_kws={'label': cbar_title})
    
    xticks = []
    xticklabels = []
    yticks = []
    yticklabels = []
    for i in range(0,MI.shape[0],args.step):
        xticks.append(i)
        xticklabels.append(i)
        yticks.append(i)
        yticklabels.append(i+args.k)
    xticks[-1] = MI.shape[0]-1; xticklabels[-1] = int(MI.shape[0]-1)
    yticks[-1] = MI.shape[0]-1+args.k; yticklabels[-1] = int(MI.shape[0]-1+args.k)
    sns_plot.set(xticks=xticks,xticklabels=xticklabels,yticks=yticks,yticklabels=yticklabels)    
    sns.despine(offset=10, trim=True)
    fig = sns_plot.get_figure()
    if args.distance=='MI': fig.savefig(args.outdir+"MI."+args.figtype,dpi=300)
    elif args.distance=='BC': fig.savefig(args.outdir+"BC."+args.figtype,dpi=300)
    elif args.distance=='BC_inv': fig.savefig(args.outdir+"BC_inv."+args.figtype,dpi=300)
    elif args.distance=='HE': fig.savefig(args.outdir+"HE."+args.figtype,dpi=300)
    elif args.distance=='JS': fig.savefig(args.outdir+"JS."+args.figtype,dpi=300)
    elif args.distance=='JS_inv': fig.savefig(args.outdir+"JS_inv."+args.figtype,dpi=300)
    end = time()
    if args.v>0: print("Plotting done in "+str(end-start)+" seconds.")
#end

plotMI()
