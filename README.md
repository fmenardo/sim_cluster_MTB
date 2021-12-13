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

```
usage: MTB_cluster_sim.py [-h] [-l lineages] [-ts time] [-br B_R] [-dr D_R] [-sr S_R] [-er E_R] [-sim_n SIM] [-cr C_R [C_R ...]] [-ps_sr [PS_SR ...]] [-ps_sy [PS_SY ...]] [-c] [-s lineages|time] [-min_mt] [-max_mt] [-t] [-rpt] [-rrt] [-SNP_t]

optional arguments:
  -h, --help            show this help message and exit
  -l lineages, --lineages lineages
                        number of infectious individuals, when the simulation exceed this number it stops
  -ts time, --time_sampling time
                        number of years of sampling (starting from present and going backward
  -br B_R, --birth_rate B_R
                        transmission rate
  -dr D_R, --death_rate D_R
                        death rate
  -sr S_R, --sampling_rate S_R
                        sampling rate
  -er E_R, --exposed_rate E_R
                        rate at which exposed become infectious
  -sim_n SIM, --simulation_number SIM
                        simulation number ID
  -cr C_R [C_R ...], --clock_rate C_R [C_R ...]
                        clock rate (nucleotide substitution per site per year), multiple values possible
  -ps_sr [PS_SR ...], --post_sim_sampling_rates [PS_SR ...]
                        probability of each strain to be sampled (post simulation), multiple rates possible at once: eg. <-ps_sr 1 0.5 0.1
  -ps_sy [PS_SY ...], --post_sim_sampling_years [PS_SY ...]
                        sample only in these years (post simulation), multiple scheme possible at once: eg. <-ps_sy 1,2,3 1,3,5> default (all)
  -c, --clean           delete all intermediate file, keep only clustering results and terminal branch lengths (default: False)
  -s lineages|time, --stop lineages|time
                        stop criterion, the MASTER simulation should stop when reaching a certain number of infectious existing lineages("lineages"; specified with -l) or after a certain time ("time",
                        specified with -t)(default = "lineages")
  -min_mt , --min_master_tips 
                        minimum number of tips in the tree output of MASTER to accept the simulation
  -max_mt , --max_master_tips 
                        max number of tips in the tree output of MASTER to accept the simulation
  -t , --time           Simulation time for MASTER, to be used with "time" as stop criterion. The simulation will stop after the amount of years specified with this option. If the lineage goes extinct
                        before the simulation is retained (Default = 10)
  -rpt , --raxml_pars_tree 
                        number of starting parsimony trees for raxml
  -rrt , --raxml_rand_tree 
                        number of starting random trees for raxml
  -SNP_t , --SNP_threshold 
                        clustering will be performed for all values in the interval 0-SNP_t (Default = 50)
```










