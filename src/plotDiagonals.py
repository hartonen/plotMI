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
    parser.add_argument("--xlim",help="Consider distances up to xlim (default=length of input sequences).",type=int,default=None)
    
    args = parser.parse_args()

    #read in the MI-matrix
    MI = np.loadtxt(args.MI)

    #compute mean of each diagonal
    mean_of_diag = []
    #compute max of each diagonal
    max_of_diag = []
    #set x-axis limit
    if args.xlim==None: xlim = int(MI.shape[0])
    elif args.xlim>int(MI.shape[0]): xlim = int(MI.shape[0])
    else: xlim = args.xlim
    
    for i in range(xlim):
        mean_of_diag.append(np.mean(np.diagonal(MI,offset=i)))
        max_of_diag.append(np.max(np.diagonal(MI,offset=i)))
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

    #plot max of each diagonal
    x = [args.k]
    for i in range(1,len(max_of_diag)): x.append(x[-1]+1)
    plt.plot(x,max_of_diag,'k')
    plt.xlim([x[0],x[-1]])
    plt.ylabel('Max MI')
    plt.xlabel('Distance between '+str(args.k)+'-mer distributions')
    plt.tight_layout()
    plt.savefig(args.outdir+"max_MI_of_diagonals."+args.figtype,dpi=150)
    plt.clf()
    np.savetxt(args.outdir+"max_MI_of_diagonals.txt",max_of_diag,delimiter='\t')
    
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
    
