#!/usr/bin/env python

import argparse
import csv
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score

def prc():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #PARAMETERS
    parser.add_argument("--class0preds",help="Full paths to files containing class 0 predictions (one file per classifier).",type=str,nargs='+')
    parser.add_argument("--class1preds",help="Full paths to files containing class 1 predictions (one file per classifier, classifier order must be same as for class 0).",type=str,nargs='+')
    parser.add_argument("--outdir",help="Full path to output directory.",type=str)
    parser.add_argument("--figtype",help="Figure type, pdf or png (default=png).",type=str,choices=['png','pdf'],default='png')
    args = parser.parse_args()

    #reading in predictions
    preds = [] #one list of predictions per classifier
    true_labels  = []
    first = True
    for f in args.class0preds:
        preds.append([])
        with open(f,'rt') as infile:
            r = csv.reader(infile,delimiter='\t')
            for row in r:
                preds[-1].append(float(row[0]))
                if first: true_labels.append(0)
            first = False    

    ind = 0        
    for f in args.class1preds:
        with open(f,'rt') as infile:
            r = csv.reader(infile,delimiter='\t')
            for row in r:
                preds[ind].append(float(row[0]))
                if ind==0: true_labels.append(1)
            ind += 1

    #calculating the precision-recall curves        
    for ind in range(0,len(preds)):
        auPRC = average_precision_score(true_labels,preds[ind])
        print("Classifier "+str(ind)+" AUprc="+str(auPRC))
        PRC = precision_recall_curve(true_labels,preds[ind])
        plt.plot(PRC[1],PRC[0],label="classifier "+str(ind)+" AUprc="+str(round(auPRC,5)))
        
    plt.xlabel('recall')
    plt.ylabel('precision')
    plt.legend()
    plt.savefig(args.outdir+'prc.'+args.figtype,dpi=150)
    plt.clf()
#end

prc()
