# Pipeline to simulate the evoluton of genome sequences of MTB under different conditions

The pipeline was tested on Ubuntu and CentOS.

## Installation with conda

```
conda env create --file MTB_sim.yml
```
This will create a conda environment with all the tools needed, with one exception. [Beast2](https://www.beast2.org/) needs to be installed manually, make sure to install the package [MASTER](http://tgvaughan.github.io/MASTER/) as well (follow the links for instuctions).

Now open the file MTB_cluster_sim.py with a text editor (this is the wrapper for the pipeline), the beginning of the file looks like this:


```
from Bio import SeqIO
from ete3 import Tree
import argparse
import os
import csv
import sys
import re
import subprocess
import random
import shutil


BEAST_PATH ="~/software/beast/bin/beast"
#SCRATCH_PATH= "/path/to/scratch"
SCRATCH_PATH= os.getcwd()

```
You need to edit the path to your beast launching script.

If you want you can specify a path to a folder where the intermediate files will be saved (handy when working on a cluster with a scratch partition). If you do so uncomment and edit the first line starting with SCRATCH_PATH, and comment the line after.

## Usage











