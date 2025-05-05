#!/bin/bash
# -----------------------------Name of the job-------------------------
#SBATCH --job-name=example
#-----------------------------Output files-----------------------------
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#-----------------------------Other information------------------------
#SBATCH --comment='Some comments'

#-----------------------------Required resources-----------------------
#SBATCH --time=1-10:52:30
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB 

#-----------------------------Environment, Operations and Job steps----
Rscript bglr_bayesian.R
