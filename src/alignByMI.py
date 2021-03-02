#!/usr/bin/env python

import argparse
import csv
import numpy as np
import gzip

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import pandas as pd
import logomaker as lm

def alignByMI():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #PARAMETERS
    parser.add_argument("--outdir",help="Full path to the output directory.",type=str)
    parser.add_argument("--seqs",help="Full path to the plain text input sequence file. Each sequence must be of same length and on its separate line.",type=str)
    parser.add_argument("--plotMI",help="Full path to plotMI output file.",type=str)
    parser.add_argument("--d",help="Distance for the interaction used in aligning (default=50bp).",type=int,default=50)
    parser.add_argument("--w",help="Visualized PWM width (default=6bp), added to both sides of the k-mer.",type=int,default=20)
    parser.add_argument("--figtype",help="png or pdf (default=png).",type=str,choices=['pdf','png'],default='png')
    parser.add_argument("--k",help="length of k-mer distributions used to calculate MI (default=3).",type=int,default=3)
    
    args = parser.parse_args()

    #read in the MI-distribution from plotMI output file
    #format of plotMI output is: ['#i','j','a','b','MI_ij(a,b)','P_ij(a,b)','P_i(a)','P_j(b)']
    #and the file is gzipped
    MI = {} #key = (i,j), value = {(a,b) : MI_ij(a,b)}
    with gzip.open(args.plotMI,'rt') as infile:
        r = csv.reader(infile,delimiter='\t')
        for row in r:
            if row[0].count('#')>0: continue
            i = int(row[0])
            j = int(row[1])
            a = row[2]
            b = row[3]
            MI_ijab = float(row[4])
            if (i,j) not in MI: MI[(i,j)] = {(a,b):MI_ijab}
            else: MI[(i,j)][(a,b)] = MI_ijab

    #read in the sequences and for each sequence select the positions with highest pointwise MI at distance d apart from each other.
    pwms = [np.zeros(shape=(4,2*args.w+args.k)),np.zeros(shape=(4,2*args.w+args.k))]
    
    count = 0.0
    with open(args.seqs,'rt') as infile:
        r = csv.reader(infile,delimiter='\t')
        for row in r:
            seq = row[0].upper()
            J = len(seq)
            seq1 = ""
            seq2 = ""
            mi_max = 0.0
            for m in range(args.w,J-2*args.k+1-args.d-args.w):
                n = m+args.d
                mi = MI[(m,n)][seq[m:m+args.k],seq[n:n+args.k]]
                if mi>mi_max:
                    seq1 = seq[m-args.w:m+args.k+args.w]
                    seq2 = seq[n-args.w:n+args.k+args.w]
                    mi_max = mi
            #Add to the PWMs the counts corresponding to the sequences coming from maximum mutual information pair
            count += 1.0
            ind = 0
            print(seq1[args.w:args.w+args.k]+"\t"+seq2[args.w:args.w+args.k])
            for s in [seq1,seq2]:
                for i in range(0,len(s)):
                    l = s[i]
                    if l=='A': pwms[ind][0,i] += 1.0
                    elif l=='C': pwms[ind][1,i] += 1.0
                    elif l=='G': pwms[ind][2,i] += 1.0
                    else: pwms[ind][3,i] += 1.0
                ind += 1 
    #visualize and save the pwms
    for ind in range(0,len(pwms)):
        plotLogo(pwms[ind],args.outdir+"pwm_"+str(ind)+"."+args.figtype)
        np.savetxt(args.outdir+"pwm_"+str(ind)+".txt",pwms[ind],delimiter='\t')
        

#end

def plotLogo(matrix,outfile):
    #plots a logo of matrix using logoMaker

    fig,ax = plt.subplots()
    
    matrix_df = pd.DataFrame(matrix.transpose())
    matrix_df.columns = ['A','C','G','T']
    logo = lm.Logo(df=matrix_df,color_scheme='classic')
    
    plt.savefig(outfile,dpi=150)
    plt.close(fig)
    plt.clf()
    plt.cla()
    #print("done!")
    return True

alignByMI()
                
                            
