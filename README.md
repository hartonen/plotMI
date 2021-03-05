# PlotMI

PlotMI is a program for visualization of pairwise interactions in a set of input sequences by computing the pairwise mutual information between positional k-mer distributions. Description of the method is in manuscript:

Hartonen T., Kivioja T., & Taipale J., (2021), "PlotMI: visualization of pairwise interactions and positional preferences learned by a deep learning model from sequence data", manuscript under preparation.


## 1. Installation

All scripts are pure Python, so no compiling is needed. Easiest way is to clone the repository to your own computer. In a desired location, type:

`git clone https://github.com/hartonen/plotMI.git`

The scripts in this repository need specific Python packages to function properly. The easiest way to make sure everything works is to create a cirtual environment [https://docs.python.org/3/library/venv.html#module-venv] containing the tested versions of each package and then run the scripts in this environment. This is done by first creating a new virtual environment:

`python3 -m venv /path/to/new/virtual/environment`

Then one needs to install all the required packages. These are listed in `data/plotMI_requirements.txt`. So activate the new virtual environment:

`source /path/to/new/virtual/environment/bin/activate`

and install the packages with pip:

`pip install -r data/plotMI_requirements.txt`

## 2. Usage examples

### 2.1 Simple example with pre-filtered input data

### 2.2 Reproducing figure 1b from plotMI-manuscript

In the following we show how to replicate the figures 1b and 1c from the plotMI-manuscript. For this we will need two other Python scripts from the authors from,  [https://github.com/hartonen/randomReads](randomReads) and  [https://github.com/hartonen/promoterAnalysis](promoterAnalysis) repositories. We will also use the  [https://bioinf.shenwei.me/seqkit/](Seqkit) tool for manipulating fasta-files.

First we generate 1 million random DNA sequences with uniform nucleotide background using a script made for this purpose:

`randomReads.py random_uniform_L200_N1M.fasta --L 200 --N 1000000`

Next we use the pre-trained CNN model to score all these sequences:

`scorePromoters.py --outfile model-36-0.990.h5-random_uniform_L200_N1M_preds.txt --model model-36-0.990.h5 --sequences random_uniform_L200_N1M.fasta --nproc 30`

Then select the fasta-IDs of sequences that score higher than 0.9 according to the model:

`awk '$2>0.9' model-36-0.990.h5-random_uniform_L200_N1M_preds.txt | cut -f1,1 > model-36-0.990.h5-random_uniform_L200_N1M_prob_more_09_ids.txt`

Then we use the wonderful seqkit-package to extract the sequences matching to these IDs and convert them to sequence-only input for mutual information plotting:

`seqkit grep -n -f model-36-0.990.h5-random_uniform_L200_N1M_prob_more_09_ids.txt random_uniform_L200_N1M.fasta | seqkit seq --seq > model-36-0.990.h5-random_uniform_L200_N1M_prob_more_09.seq`

These sequences can then be visualized using plotMI (Figure 1b):

`plotMI.py --outdir model-36-0.990.h5-random_uniform_L200_N1M_prob_more_09- --seqs model-36-0.990.h5-random_uniform_L200_N1M_prob_more_09.seq --nproc 30 --figtype png --k 3 --v 1 --p 5`

Note that due to generating the input sequences by random, the figure will not look exactly the same as in the manuscript.
