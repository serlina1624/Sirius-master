#!/bin/bash

#SBATCH --job-name=trinity
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=25
#SBATCH -p normal
#SBATCH -t 11520 
#SBATCH -o results/LOGS/trinity.out 
#SBATCH -e results/LOGS/trinity.err 
#SBATCH -c 12 

# Your script goes here
source ~/.bash_profile
/home/telnaya/programm/TRINITY/Trinity --seqType fa --single /home/telnaya/data/All_spicies.fasta --max_memory 200G  --CPU 25 --output /home/telnaya/results/trinity --verbose


