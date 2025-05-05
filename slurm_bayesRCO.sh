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

bayesRCO_path/bayesRCO -bfile mice_835_qc_GRCm39_ref -out mice_835_qc_GRCm39_ref -ncat 1 -catfile annot.txt -burnin 20000 -numit 120000 -thin 100

bayesRCO_path/bayesRCO -bfile mice_835_qc_GRCm39_valid -predict -out mice_835_qc_GRCm39_valid -model mice_835_qc_GRCm39_ref.model -freq mice_835_qc_GRCm39_ref.frq -param mice_835_qc_GRCm39_ref.param -ncat 1 -catfile annot.txt