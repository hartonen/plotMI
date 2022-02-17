#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import seaborn as sns

import argparse
import csv
import numpy as np
import gzip
import sys

from statsmodels.stats.multitest import multipletests

def calcPvalues():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #PARAMETERS
    parser.add_argument("--outdir",help="Full path to the output directory.",type=str)
    parser.add_argument("--infiles",help="Full path to a plain text input file. Each line contains path to a single file that contains an MI matrix, as output by plotMI. First line is assumed to contain the file with MI matrix from sequences filtered by model score. Other lines are assumed to contain files with MI matrices computed from random samples.",type=str)
    parser.add_argument("--nproc",help="Number of parallel processes used when computing MI (default=1).",type=int,default=1)
    parser.add_argument("--figtype",help="png or pdf (default=png).",type=str,choices=['pdf','png'],default='png')
    parser.add_argument("--k",help="length of k-mer distributions used to calculate MI (default=3).",type=int,default=3)
    parser.add_argument("--v",help="Verbosity level, 0=none, 1=print info on screen (default=1).",type=int,choices=[0,1],default=1)
    parser.add_argument("--p",help="Significance threshold (default=0.05).",type=float,default=0.05)
    parser.add_argument("--colorscale",help="If set to log, colormap is scaled logarithmically (default=lin, meaning linear scaling).",type=str,choices=['lin','log'],default='lin')
    parser.add_argument("--minmi",help="Set minimum value for colormap, helpful if you want to be sure that the minimum value is 0 (default=minimum value in MI matrix).",default=None,type=float)
    parser.add_argument("--step",help="Step size for axis ticks in MI-plot (default=20).",type=int,default=20)
    parser.add_argument("--method",help="Method for adjusting p-values for multiple hypothesis testing. Should be one accepted by statsmodels.stats.multitest (see: https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html). Default is fdr_bh for Benjamini-Hochberg procedure.",type=str,default='fdr_bh')

    args = parser.parse_args()

    #save the command used to evoke calcPvalues to a log file
    with open(args.outdir+"calcPvalues_log.txt",'wt') as logfile:
        logfile.write(" ".join(sys.argv))

    #read in the file names
    filenames = []
    with open(args.infiles,'rt') as infile:
        for row in infile: filenames.append(row.strip())
        
    #read in the MI computed from model-filtered sequences
    MI_model = np.loadtxt(filenames[0])
    #read in the MIs computed from random samples
    MI_random = np.zeros(shape=(MI_model.shape[0],MI_model.shape[1],len(filenames)-1))
    for i in range(1,len(filenames)): MI_random[:,:,i-1] = np.loadtxt(filenames[i])

    #for each position in MI_model, check how many random samples have a higher MI value
    Ps = np.zeros(shape=MI_model.shape)
    for i in range(0,len(filenames)-1): Ps[np.where(MI_random[:,:,i]>MI_model)] += 1.0
    Ps /= (len(filenames)-1)

    #plot the p-values
    np.savetxt(args.outdir+"p-vals.txt.gz",Ps,delimiter='\t')
    plotMatrix(Ps,'empirical p-value',args.outdir+"p-vals."+args.figtype,args,colorscale='log')

    #adjust for multiple hypothesis testing using any method available in statsmodels
    iu = np.triu_indices(Ps.shape[0])
    multi = multipletests(Ps[iu],alpha=args.p,method=args.method)[1]
    Ps_adj = np.zeros(shape=Ps.shape)
    Ps_adj[iu] = multi
    #plot significant vs non-significant pairs
    Ps_adj = Ps_adj+Ps_adj.T-np.diag(np.diag(Ps_adj)) #copy upper triangle values to lower triangle
    significant = np.zeros(shape=Ps.shape)
    significant[np.where(Ps_adj<=args.p)] = 1.0
    np.savetxt(args.outdir+"significant_vs_non_significant.txt.gz",significant,delimiter='\t')
    plotMatrix(significant,'significant MI signal in black (p_adj<0.05)',args.outdir+"significant_vs_non_significant."+args.figtype,args,cmap='Greys',binaryticks=True)

    
    #plot MI_model-np.mean(MI_random)
    MI_model_minus = MI_model-np.mean(MI_random,axis=2)
    np.savetxt(args.outdir+"MI_model-MI_mean_from_random_samples.txt.gz",MI_model_minus,delimiter='\t')
    plotMatrix(MI_model_minus,'MI_model-<MI_random>',args.outdir+"MI_model-MI_mean_from_random_samples."+args.figtype,args)

#end

def plotMatrix(matrix,cbar_title,outname,args,cmap='viridis',colorscale='linear',colorbar=True,binaryticks=False):
    #matrix = matrix to be plotted
    #cbar_title = colorbar title
    #outname = full path to output file
    #args = command line arguments
    
    if args.minmi!=None:
        if colorscale=='log': sns_plot = sns.heatmap(matrix,cmap=cmap,cbar=colorbar,cbar_kws={'label': cbar_title},vmin=args.minmi,norm=LogNorm())
        else: sns_plot = sns.heatmap(matrix,cmap=cmap,cbar=colorbar,cbar_kws={'label': cbar_title},vmin=args.minmi)
    else:
        if colorscale=='log': sns_plot = sns.heatmap(matrix,cmap=cmap,cbar=colorbar,cbar_kws={'label': cbar_title},norm=LogNorm())
        else: sns_plot = sns.heatmap(matrix,cmap=cmap,cbar=colorbar,cbar_kws={'label': cbar_title})

    xticks = []
    xticklabels = []
    yticks = []
    yticklabels = []
    for i in range(0,matrix.shape[0],args.step):
        xticks.append(i)
        xticklabels.append(i)
        yticks.append(i)
        yticklabels.append(i+args.k)
    xticks[-1] = matrix.shape[0]-1; xticklabels[-1] = int(matrix.shape[0]-1)
    yticks[-1] = matrix.shape[0]-1+args.k; yticklabels[-1] = int(matrix.shape[0]-1+args.k)
    sns_plot.set(xticks=xticks,xticklabels=xticklabels,yticks=yticks,yticklabels=yticklabels)
    sns.despine(offset=10, trim=True)

    if binaryticks:
        cbar = sns_plot.collections[0].colorbar
        cbar.set_ticks([0,1])
        cbar.set_ticklabels(['p>='+str(args.p), 'p<'+str(args.p)])
    
    fig = sns_plot.get_figure()
    fig.savefig(outname,dpi=300)
    plt.clf()
    
calcPvalues()
