#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sns

from time import time
import argparse
import csv
import numpy as np

import multiprocessing as mp

from helpers import getP_j, getMI_mn

def plotMI():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #PARAMETERS
    parser.add_argument("--outdir",help="Full path to the output directory.",type=str)
    parser.add_argument("--seqs",help="Full path to the plain text input sequence file. Each sequence must be of same length and on its separate line.",type=str)
    parser.add_argument("--nproc",help="Number of parallel processes used when computing MI (default=1).",type=int,default=1)
    parser.add_argument("--figtype",help="png or pdf (default=png).",type=str,choices=['pdf','png'],default='png')
    parser.add_argument("--k",help="length of k-mer distributions used to calculate MI (default=3).",type=int,default=3)
    parser.add_argument("--v",help="Verbosity level, 0=none, 1=print info on screen (default=1).",type=int,choices=[0,1],default=1)
    parser.add_argument("--p",help="Pseudocount added to each k-mer count (default=0).",type=int,default=0)
    args = parser.parse_args()
    
    #read in the sequences and store them as strings
    start = time()
    seqs = []
    with open(args.seqs,'rt') as infile:
        r = csv.reader(infile)
        for row in r:
            if len(row)<1: continue
            seqs.append(row[0])
    
    I = len(seqs) #number of sequences
    J = len(seqs[0]) #length of sequences
    
    end = time()
    if args.v>0: print("Read in "+str(I)+" sequences of length "+str(J)+" in "+str(end-start)+" seconds.")
    
    start = time()
    #calculate single site frequencies of k-mers in parallel
    pool = mp.Pool(args.nproc)
    res = [pool.apply_async(getP_j,args=(seqs,j,I,J,args.k,args.p)) for j in range(0,J-args.k+1)]
    P = [r.get() for r in res] #list containing the single site k-mer frequencies. Each distribution is saved as a dictionary where k-mer is key and its frequency is value
    
    #print(P[0])
    #print(sum([P[0][kmer] for kmer in P[0]]))
    pool.close()
    pool.join()
    pool.terminate()
    
    end = time()
    if args.v>0: print("Calculated the positional "+str(args.k)+"-mer frequencies in "+str(end-start)+" seconds.")
    
    start = time()
    #calculate MI in parallel
    pool = mp.Pool(args.nproc)
    res = []
    
    M = 0 #number of pairwise k-mer distributions
    for m in range(0,J-2*args.k+1):
        M += 1
        #overlapping_inds = set([l for l in range(m-args.k,m+args.k)])
        for n in range(m+args.k,J-args.k+1):
            #if n in overlapping_inds: continue
            #print("m="+str(m)+", n="+str(n))
            res.append(pool.apply_async(getMI_mn,args=(seqs,m,n,I,J,P,args.k,args.p)))
    res = [r.get() for r in res]    

    pool.close()
    pool.join()
    pool.terminate()
    
    #build the MI matrix
    #print("M="+str(M))
    MI = np.zeros(shape=(M,M))
    #print(MI.shape)
    for j in range(0,len(res)):
        inds = res[j][0]
        #print(inds)
        mi = res[j][1]
        MI[inds[1]-args.k,inds[0]] = mi

    xticks = [res[0][0][1],res[-1][0][1]]
    yticks = [res[0][0][0],res[-1][0][0]]
    end = time()
    if args.v>0: print("Computed mutual information in "+str(end-start)+" seconds.")
    
    start = time()
    #save the MI matrix to file
    np.savetxt(args.outdir+"MI.txt.gz",MI,delimiter='\t')
    
    #plot the lower triangle of the MI-matrix
    mask = np.zeros_like(MI)
    mask[np.triu_indices_from(mask)] = True
    print(MI)
    sns_plot = sns.heatmap(MI,cmap='Blues',cbar=True,cbar_kws={'label': 'MI'},mask=mask,xticklabels=xticks,yticklabels=yticks)
    sns_plot.set_xticks([i-args.k for i in xticks])
    sns_plot.set_yticks(yticks)
    fig = sns_plot.get_figure()
    #fig.set_xlabel('position')
    #fig.set_ylabel('position')
    fig.savefig(args.outdir+"MI."+args.figtype,dpi=150)
    
    end = time()
    if args.v>0: print("Plotting done in "+str(end-start)+" seconds.")
#end

plotMI()
