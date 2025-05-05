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

bayesR_path/bayesR -bfile mice_835_qc_GRCm39 -out mice_835_qc_GRCm39 -numit 120000 -burnin 20000 -thin 100 -n 1

bayesR_path/bayesR -bfile mice_835_qc_GRCm39 -out predict_ebv -predict -model mice_835_qc_GRCm39.model -freq mice_835_qc_GRCm39.frq -param mice_835_qc_GRCm39.param -n 1
