# PlotMI

PlotMI is a program for visualization of pairwise interactions in a set of input sequences by computing the pairwise mutual information between positional k-mer distributions. Description of the method is in manuscript:



## 1. Installation

All scripts are pure Python, so no compiling is needed. Easiest way is to clone the repository to your own computer. In a desired location, type:

`git clone https://github.com/hartonen/plotMI.git`

The scripts in this repository need specific Python packages to function properly. The easiest way to make sure everything works is to (create a virtual environment) [https://docs.python.org/3/library/venv.html#module-venv] containing the tested versions of each package and then run the scripts in this environment. This is done by first creating a new virtual environment:

`python3 -m venv /path/to/new/virtual/environment`

Then one needs to install all the required packages. These are listed in `data/plotMI_requirements.txt`. So activate the new virtual environment:

`source /path/to/new/virtual/environment/bin/activate`

and install the packages with pip:

`pip install -r data/plotMI_requirements.txt`

## 2. Usage example

Sequences in file `data/` can be used to re-create Figure 1b of the plotMI manuscript by typing: 