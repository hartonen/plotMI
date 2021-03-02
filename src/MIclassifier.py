#!/usr/bin/env python

import argparse
import csv
import numpy as np
import gzip

from math import log
from math import exp

def MIclassifier():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #PARAMETERS
    parser.add_argument("--outdir",help="Full path to the output directory.",type=str)
    parser.add_argument("--seqs",help="Full path to the plain text input sequence file for scoring. Each sequence must be of same length and on its separate line.",type=str)
    parser.add_argument("--plotMI",help="Full path to plotMI output file.",type=str)
    parser.add_argument("--k",help="length of k-mer distributions used to calculate MI (default=3).",type=int,default=3)
    args = parser.parse_args()

    #read in the MI-distribution from plotMI output file
    #format of plotMI output is: ['#i','j','a','b','MI_ij(a,b)','P_ij(a,b)','P_i(a)','P_j(b)']
    #and the file is gzipped
    MI = {} #key = (i,j), value = {(a,b) : P_ij(a,b)} 
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

    #then scoring the input sequences using the sum of MI as classifier
    with open(args.outdir+"MI_classifier_scores.txt",'wt') as outfile:
        w = csv.writer(outfile)
        with open(args.seqs,'rt') as infile:
            r = csv.reader(infile,delimiter='\t')
            for row in r:
                seq = row[0].upper()
                J = len(seq)
                score = 0.0
                for m in range(0,J-2*args.k+1):
                    for n in range(m+args.k,J-args.k+1): score += MI[m,n][seq[m:m+args.k],seq[n:n+args.k]]
                w.writerow([score])    

#end

MIclassifier()    
