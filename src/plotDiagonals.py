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
    parser.add_argument("--figtype",help="Figure type, pdf (=default) or png.",type=str,choices=['pdf','png'],default='pdf')
    
    args = parser.parse_args()

    #read in the MI-matrix
    MI = np.loadtxt(args.MI)

    #compute mean of each diagonal
    mean_of_diag = []
    for i in range(int(MI.shape[0]/2)+1): mean_of_diag.append(np.mean(np.diagonal(MI,offset=i)))
    #plot mean of each diagonal
    plt.plot(mean_of_diag,'k')
    plt.ylabel('Mean MI of diagonal')
    plt.xlabel('Offset from main diagonal')
    plt.tight_layout()
    plt.savefig(args.outdir+"mean_MI_of_diagonals."+args.figtype,dpi=150)
    plt.clf()
    np.savetxt(args.outdir+"mean_MI_of_diagonals.txt",mean_of_diag,delimiter='\t')

    #plot MI of main diagonal
    plt.plot(np.diagonal(MI,offset=0),'k')
    plt.xlabel('Position')
    plt.ylabel('MI')
    plt.title('Main diagonal')
    plt.tight_layout()
    plt.savefig(args.outdir+'main_diagonal_MI.'+args.figtype,dpi=150)
    plt.clf()
    np.savetxt(args.outdir+'main_diagonal_MI.txt',np.diagonal(MI,offset=0),delimiter='\t')

#end

plotDiagonals()    
    
