# PlotMI

PlotMI is a program for visualization of pairwise interactions in a set of input sequences by computing the pairwise mutual information between positional k-mer distributions. Description of the method is in manuscript:

Hartonen T., Kivioja T., & Taipale J., (2021), "PlotMI: visualization of pairwise interactions and positional preferences learned by a deep learning model from sequence data", manuscript under preparation.


## 1. Installation

All scripts are pure Python, so no compiling is needed. Easiest way is to clone the repository to your own computer. In a desired location, type:

`git clone https://github.com/hartonen/plotMI.git`

The scripts in this repository need specific Python packages to function properly. The easiest way to make sure everything works is to create a virtual environment (https://docs.python.org/3/library/venv.html#module-venv) containing the tested versions of each package and then run the scripts in this environment. This is done by first creating a new virtual environment:

`python3 -m venv /path/to/new/virtual/environment`

Then one needs to install all the required packages. These are listed in `data/plotMI_requirements.txt`. So activate the new virtual environment:

`source /path/to/new/virtual/environment/bin/activate`

and install the packages with pip:

`pip install -r data/plotMI_requirements.txt`

## 2. Usage

Help message can be evoked by typing:

```
plotMI.py -h
usage: plotMI.py [-h] [--outdir OUTDIR] [--seqs SEQS] [--nproc NPROC]
                 [--figtype {pdf,png}] [--k K] [--v {0,1}] [--p P]
                 [--alphabet ALPHABET] [--minmi MINMI] [--step STEP]
                 [--save_distributions {yes,no}]

optional arguments:
  -h, --help            show this help message and exit
  --outdir OUTDIR       Full path to the output directory.
  --seqs SEQS           Full path to the plain text input sequence file. Each
                        sequence must be of same length and on its separate
                        line.
  --nproc NPROC         Number of parallel processes used when computing MI
                        (default=1).
  --figtype {pdf,png}   png or pdf (default=png).
  --k K                 length of k-mer distributions used to calculate MI
                        (default=3).
  --v {0,1}             Verbosity level, 0=none, 1=print info on screen
                        (default=1).
  --p P                 Multiplier for pseudocount mass added to k-mer count.
                        Total pseudocount mass added is p*(number of
                        sequences) (default=5).
  --alphabet ALPHABET   A string containing each individual letter in the
                        alphabet used (default=ACGT). NOTE! This is case-
                        sensitive.
  --minmi MINMI         Set minimum value for colormap, helpful if you want to
                        be sure that the minimum value is 0 (default=minimum
                        value in MI matrix).
  --step STEP           Step size for axis ticks in MI-plot (default=20).
  --save_distributions {yes,no}
                        If yes, save the positional and pairwise k-mer
                        distributions and the MI contributions from each k-mer
                        and position pair into a separate file (default=no).
                        Note that this is a large file of approximately 1GB.
```

## 3. Usage examples

### 3.1 Simple example with pre-filtered input data

The most basic use case for plotMI only requires input sequences that have been filtered using a machine learning model. Here as an example, we reproduce Figure 1d from the plotMI manuscript using random uniform DNA sequences that have been filtered using a convolutional neural network model that has been trained to recognize human promoter sequences (see the manuscript for details). The file `data/model-492-0.882.h5-random-uniform-1M-prob-more-09.seq` contains 49,261 sequences that have a probability >0.9 of being active human promoters based on the CNN model. Using these sequences, we can visualize the interactions learned by the human promoter CNN model with:  

`plotMI.py --outdir ./test- --seqs model-492-0.882.h5-random-uniform-1M-prob-more-09.seq --nproc 4 --figtype png
Alphabet: ACTG
Read in 49261 sequences of length 100 in 0.04489850997924805 seconds.
Calculated the positional 3-mer frequencies in 2.0755910873413086 seconds.
Computed mutual information in 83.16843390464783 seconds.
Plotting done in 0.4934566020965576 seconds.
`

This should produce two output files identical to files `data/test-MI.png` and `data/test-MI.txt.gz`.

### 3.2 Reproducing figure 1b from plotMI-manuscript

In the following we show how to replicate the figures 1b and 1c from the plotMI-manuscript. For this we will need two other Python scripts from the authors from,  [https://github.com/hartonen/randomReads](randomReads) and  [https://github.com/hartonen/promoterAnalysis](promoterAnalysis) repositories. We will also use the  [https://bioinf.shenwei.me/seqkit/](Seqkit) tool for manipulating fasta-files.

First we generate 10 million random DNA sequences (10 separate files for faster scoring with the CNN model) with uniform nucleotide background using a script made for this purpose:

`for i in {1..10}; do randomReads.py random_uniform_L200_N1M_sample"$i".fasta --L 200 --N 1000000; done;`

Next we use the pre-trained CNN model to score all these sequences. The pre-trained model can be downloaded from Zenodo with DOI: 10.5281/zenodo.4596516.

`for i in {1..10}; do scorePromoters.py --outfile model-36-0.990.h5-random_uniform_L200_N1M_sample"$i"_preds.txt --model model-36-0.990.h5 --sequences random_uniform_L200_N1M_sample"$i".fasta --nproc 4 & done;`

Then select the fasta-IDs of sequences that score higher than 0.9 according to the model:

`for i in {1..10}; do awk '$2>0.9' model-36-0.990.h5-random_uniform_L200_N1M_sample"$i"_preds.txt | cut -f1,1 > model-36-0.990.h5-random_uniform_L200_N1M_sample"$i"_prob_more_09_ids.txt; done;`

Then we use the wonderful seqkit-package to extract the sequences matching to these IDs and convert them to sequence-only input for mutual information plotting:

`for i in {1..10}; do seqkit grep -n -f model-36-0.990.h5-random_uniform_L200_N1M_sample"$i"_prob_more_09_ids.txt random_uniform_L200_N1M_sample"$i".fasta | seqkit seq --seq > model-36-0.990.h5-random_uniform_L200_N1M_sample"$i"_prob_more_09.seq; done;`

Combine the sequences to a single file for MI plotting:

`cat model-36-0.990.h5-random_uniform_L200_N1M_sample*_prob_more_09.seq > model-36-0.990.h5-random_uniform_L200_N10M_prob_more_09.seq`

These sequences can then be visualized using plotMI (Figure 1b):

`plotMI.py --outdir model-36-0.990.h5-random_uniform_L200_N1M0_prob_more_09- --seqs model-36-0.990.h5-random_uniform_L200_N10M_prob_more_09.seq --nproc 16 --figtype png --k 3 --v 1 --p 5`

Note that due to generating the input sequences by random, the figure will not look exactly the same as in the manuscript.

## 4. Output description

By default plotMI will output two files: an image that contains the MI plot and a gzipped tab-delimited text-file that contains the corresponding MI matrix. Optionally, one can set the flag `--save_distribution yes`, which will then output the estimated probabilities and MI contributions for each k-mer and position pair. Note that depending on the length of the model, this is a large file (>1GB) An example of ten first lines of this file:

```
#i      j       a       b       MI_ij(a,b)      P_ij(a,b)       P_i(a)  P_j(b)
0       3       GGG     GGG     -3.475747040270734e-05  0.00020345052083333334  0.01510622461859291     0.015161103336626056
0       3       GGG     GGC     4.352340047630113e-05   0.0002583292388664801   0.01510622461859291     0.015215982054659204
0       3       GGG     GGT     -4.621769658297288e-05  0.00020345052083333334  0.01510622461859291     0.01576476923499067
0       3       GGG     GGA     -4.926708687151659e-05  0.00020345052083333334  0.01510622461859291     0.01592940538909011
0       3       GGG     GCG     0.00012695982387922662  0.00031320795689962684  0.01510622461859291     0.015655011798924378
0       3       GGG     GCC     2.516479878216178e-05   0.0002583292388664801   0.01510622461859291     0.015984284107123256
0       3       GGG     GCT     -3.155270299332963e-05  0.00020345052083333334  0.01510622461859291     0.014996467182526617
0       3       GGG     GCA     2.516479878216178e-05   0.0002583292388664801   0.01510622461859291     0.015984284107123256
0       3       GGG     GTG     4.893921467700538e-05   0.0002583292388664801   0.01510622461859291     0.014996467182526617

```
