#!/bin/bash
# -----------------------------Name of the job-------------------------
#SBATCH --job-name=example_asreml
#-----------------------------Output files-----------------------------
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#-----------------------------Other information------------------------
#SBATCH --comment='Some comments'

#-----------------------------Required resources-----------------------
#SBATCH --time=1-10:52:30
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5GB

#-----------------------------Environment, Operations and Job steps----
module load WUR/ABGC/asreml
asreml asreml_gblup.as