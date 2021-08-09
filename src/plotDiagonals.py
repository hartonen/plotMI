#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import argparse
import numpy as np

def plotDiagonals():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #PARAMETERS
    parser.add_argument("--outdir",help="Full path to the output directory.",type=str)
    parser.add_argument("--MI",help="Full path to the tab-separated txt-file containing the MI-matrix.",type=str)
    parser.add_argument("--k",help="Value of k used to create the MI-matrix.",type=int,default=3)
    parser.add_argument("--figtype",help="Figure type, pdf (=default) or png.",type=str,choices=['pdf','png'],default='pdf')
    
    args = parser.parse_args()

    #read in the MI-matrix
    MI = np.loadtxt(args.MI)

    #compute mean of each diagonal
    mean_of_diag = []
    for i in range(int(MI.shape[0]/2)+1): mean_of_diag.append(np.mean(np.diagonal(MI,offset=i)))
    #plot mean of each diagonal
    x = [args.k]
    for i in range(1,len(mean_of_diag)): x.append(x[-1]+1)
    plt.plot(x,mean_of_diag,'k')
    plt.xlim([x[0],x[-1]])
    plt.ylabel('Mean MI')
    plt.xlabel('Distance between '+str(args.k)+'-mer distributions')
    plt.tight_layout()
    plt.savefig(args.outdir+"mean_MI_of_diagonals."+args.figtype,dpi=150)
    plt.clf()
    np.savetxt(args.outdir+"mean_MI_of_diagonals.txt",mean_of_diag,delimiter='\t')

    #plot MI of main diagonal
    main_diag = np.diagonal(MI,offset=0)
    x = [args.k]
    for i in range(1,len(main_diag)): x.append(x[-1]+1)
    plt.plot(main_diag,'k')
    plt.xlabel('Position')
    plt.ylabel('MI')
    #plt.title('Main diagonal')
    plt.tight_layout()
    plt.savefig(args.outdir+'main_diagonal_MI.'+args.figtype,dpi=150)
    plt.clf()
    np.savetxt(args.outdir+'main_diagonal_MI.txt',np.diagonal(MI,offset=0),delimiter='\t')

#end

plotDiagonals()    
    
