# Pipeline to simulate the evolution of genome sequences of MTB under different conditions

This repository containes supplementary data and code used for the preprint **"Clustering and terminal branch lengths analyses can not reliably detect variation in transmission among sub-populations of *Mycobacterium tuberculosis*"**. 


* sim_results:        a folder containing results for all simulations performed in the study (clustering rates and terminal branch lengths)
* MTB_cluster_sim.py: the python wrapper fro the pipeline. The pipeline was tested on Ubuntu and CentOS.
* collect_res.R:      R script to join outputs (after running several simulations with the same settings)
* plot_results.R:     R script to compare different settings, produces plots and summary table 
* MTB_sim.yml:        yml file to create the conda environment to run the pipeline
* R_info.yml:         yml file yml file to create the conda environment to run the R scripts


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
usage: MTB_cluster_sim.py [-h] [-l lineages] [-ts time] [-br B_R] [-dr D_R] [-sr S_R] [-er E_R] [-sim_n SIM] 
        [-cr C_R [C_R ...]] [-ps_sr [PS_SR ...]] [-ps_sy [PS_SY ...]] [-c] [-s lineages|time]
        [-min_mt] [-max_mt] [-t] [-rpt] [-rrt] [-SNP_t]

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
                        probability of each strain to be sampled (post simulation), 
                        multiple rates possible at once: eg. <-ps_sr 1 0.5 0.1
  -ps_sy [PS_SY ...], --post_sim_sampling_years [PS_SY ...]
                        sample only in these years (post simulation), multiple scheme possible at once: 
                        eg. <-ps_sy 1,2,3 1,3,5> default (all)
  -c, --clean           delete all intermediate file, keep only clustering results and terminal branch lengths 
                        (default: False)
  -s lineages|time, --stop lineages|time
                        stop criterion, the MASTER simulation should stop when reaching a certain number 
                        of infectious existing lineages("lineages"; specified with -l) or after a certain time 
                        ("time" specified with -t)(default = "lineages")
  -min_mt , --min_master_tips 
                        minimum number of tips in the tree output of MASTER to accept the simulation
  -max_mt , --max_master_tips 
                        max number of tips in the tree output of MASTER to accept the simulation
  -t , --time           Simulation time for MASTER, to be used with "time" as stop criterion. 
          The simulation will stop after the amount of years specified with this option. If the lineage goes extinct
                        before the simulation is retained (Default = 10)
  -rpt , --raxml_pars_tree 
                        number of starting parsimony trees for raxml
  -rrt , --raxml_rand_tree 
                        number of starting random trees for raxml
  -SNP_t , --SNP_threshold 
                        clustering will be performed for all values in the interval 0-SNP_t (Default = 50)
```

## Stop conditions

There are two possible stop conditions which are mutually exclusive: `-s lineages|time`. 

`lineages`: the simulation stops when a number of living infectious individuals is reached (specified with `-l`). If the lineage goes extinct before the simulation is discarded. 

`time`: the simulation stop after a certain number of years (specified with `-t`). If the lineage goes extinct before the simulation is retained.

## An example

If we want to simulate the same conditions corresponding to sup-population type 2 in the manuscript (Fig. 2)

```
conda activate MTB_sim_env
python MTB_cluster_sim.py -ts 10 -br 1.53 -dr 0.85 -sr 0.85 -er 1.7 -sim_n 1 -cr 0.00000007 -s time -t 30 -min_mt 100 --clean
```

## Output

Output files (these files are not deleted with `--clean`):

* .tr.newick:      the transmission tree simulated by MASTER
* .newick:         the corresponding phylogenetic tree
* .fasta_var:      the fasta file with the simulated sequences (only variable positions)
* .raxml_rescaled: the inferred phylogenetic tree (with corrected branch length)
* .cl_r:           the clustering rates at different SNP thresholds (from 0 to `SNP_t`)
* _ldist.csv:       the list of terminal branch lengths in SNPs (tip_ID,SNP)

The name of the output files (before extension) is composed by the values of the different parameters:

`br`\_`er`\_`dr`\_`sr`\_`l`\_`min_mt`-`max_mt`\_`ts`\_`t`\_`cr`\_`er`\_`ps_sy`\_`ps_sr`\_`sim_n`

## Collect results

There is a R script available for post-processing of simulation results, it relies on `R.utils`.

I have a conda environment for this and for the plotting script:

You might need to run this:

```conda config --append channels conda-forge```
and then 
```conda create --name R_info --file R_info.yml```


If you ran multiple simulations with the same settings eg.(will take a few minutes, decrease number of simulations or `min_mt` for testing ):

```
for i in {1..10}; do python MTB_cluster_sim.py -ts 10 -br 1.53 -dr 0.85 -sr 0.85 -er 1.7 -sim_n $i -cr 0.00000007 -s time -t 30 -min_mt 100 --clean ;echo "this is the $i simulation"; done;
```
You can collect the results of the 10 simulations as follow:

```
conda activate R_info
cd sim_1.53_1.7_0.85_0.85_0_100-2500_10_30/
Rscript ../collect_res.R 1.53_1.7_0.85_0.85_0_100-2500_10_30_7e-08_1,2,3,4,5,6,7,8,9,10_1
```

This will generate 4 files:

* `1.53_1.7_0.85_0.85_0_100-2500_10_30_7e-08_1,2,3,4,5,6,7,8,9,10_1.all_cl_r_concat`      : clustering rates for all simulations (one sim per line)
* `1.53_1.7_0.85_0.85_0_100-2500_10_30_7e-08_1,2,3,4,5,6,7,8,9,10_1.all_ldist_concat.csv` : terminal branch lengths for all tips in all simulations (one tip per line)
* `1.53_1.7_0.85_0.85_0_100-2500_10_30_7e-08_1,2,3,4,5,6,7,8,9,10_1.info`                 : some general info and stats
* `1.53_1.7_0.85_0.85_0_100-2500_10_30_7e-08_1,2,3,4,5,6,7,8,9,10_1.N_count`              : the number of tips in each simulated dataset (one sim per line)

## Plot results

There is a R script available for plotting the results of different simulation settings, it uses the following packages: `argparser`, `ggplot2`, `data.table`, `ggpubr`.

here the usage:

```
usage: plot_results.R [--] [--help] [--opts OPTS] [-S -S [-f -F [-l -L
       [-o -O [-i -I

compare different simulation settings, generate plots and a summary res
table

flags:
  -h, --help  show this help message and exit

optional arguments:
  -x, --opts  RDS file containing argument values
  -S, -S      plot clustering rates for all SNP threshold =< to this
              value [default: 5]
  -f, -f      stem for output files
  -l, -l      x labels for plots. If spaces in label use quotes(eg. -l
              "x label")
  -o, -o      comma delimited list of labels for each file in the
              desired order, eg. -o 1,2,3,4. If there are spaces within
              label use quotes
  -i, -i      colon delimited list of stem name for each simulation
              (including path) same order as in -o. Use quotes if
              spaces in path

```



To reproduce Figure 1 (will take a few minutes):

```
conda activate R_info
cd sim_results
Rscript ../plot_results.R -S 7 -f Figure1 -l "Median infectious period (months) - R0" -o "17 - 0.9,17 - 1,17 - 1.1,8 - 0.9,8 - 1,8 - 1.1,4 - 0.9,4 - 1,4 - 1.1" -i 0.45_1.0_0.0_0.5_0_100-2500_10_30_8e-08_1,2,3,4,5,6,7,8,9,10_1:0.5_1.0_0.0_0.5_0_100-2500_10_30_8e-08_1,2,3,4,5,6,7,8,9,10_1:0.55_1.0_0.0_0.5_0_100-2500_10_30_8e-08_1,2,3,4,5,6,7,8,9,10_1:0.9_1.0_0.0_1.0_0_100-2500_10_30_8e-08_1,2,3,4,5,6,7,8,9,10_1:1.0_1.0_0.0_1.0_0_100-2500_10_30_8e-08_1,2,3,4,5,6,7,8,9,10_1:1.1_1.0_0.0_1.0_0_100-2500_10_30_8e-08_1,2,3,4,5,6,7,8,9,10_1:1.8_1.0_0.0_2.0_0_100-2500_10_30_8e-08_1,2,3,4,5,6,7,8,9,10_1:2.0_1.0_0.0_2.0_0_100-2500_10_30_8e-08_1,2,3,4,5,6,7,8,9,10_1:2.2_1.0_0.0_2.0_0_100-2500_10_30_8e-08_1,2,3,4,5,6,7,8,9,10_1
```

To reproduce Figure 2:

```
conda activate R_info
cd sim_results
Rscript ../plot_results.R -S 9 -f Figure2 -l "" -o "Type 1,Type 2" -i 0.77_0.7_0.35_0.35_0_100-2500_10_30_1e-07_1,2,3,4,5,6,7,8,9,10_1:1.53_1.7_0.85_0.85_0_100-2500_10_30_7e-08_1,2,3,4,5,6,7,8,9,10_1

```

