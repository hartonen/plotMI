#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sns

from time import time
import argparse
import csv
import numpy as np
from numpy import ma

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
    parser.add_argument("--p",help="Multiplier for pseudocount mass added to k-mer count. Total pseudocount mass added is p*(number of sequences) (default=5).",type=int,default=5)
    parser.add_argument("--alphabet",help="A string containing each individual letter in the alphabet used (default=ACGT). NOTE! This is case-sensitive.",type=str,default="ACGT")
    args = parser.parse_args()
    
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
    #calculate MI in parallel
    pool = mp.Pool(args.nproc)
    res = []
    
    M = 0 #number of pairwise k-mer distributions
    for m in range(0,J-2*args.k+1):
        M += 1
        for n in range(m+args.k,J-args.k+1): res.append(pool.apply_async(getMI_mn,args=(seqs,m,n,I,J,P,args.k,p,alphabet)))
    res = [r.get() for r in res]    

    pool.close()
    pool.join()
    pool.terminate()
    
    #build the MI matrix
    MI = -1*np.ones(shape=(M,M))
    for j in range(0,len(res)):
        inds = res[j][0]
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
    mask[MI<0] = True

    sns_plot = sns.heatmap(MI,cmap='jet',cbar=True,cbar_kws={'label': 'MI'},xticklabels=xticks,yticklabels=yticks,mask=mask,vmin=0)
    sns_plot.set_xticks([i-args.k for i in xticks])
    sns_plot.set_yticks(yticks)
    sns.despine(offset=10, trim=True)
    fig = sns_plot.get_figure()
    fig.savefig(args.outdir+"MI."+args.figtype,dpi=150)
    
    end = time()
    if args.v>0: print("Plotting done in "+str(end-start)+" seconds.")
#end

plotMI()
