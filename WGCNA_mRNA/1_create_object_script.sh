#!/bin/bash
#SBATCH --job-name=wgcna
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p g100_usr_prod
#SBATCH --mem=128GB
#SBATCH --time 5:59:00
#SBATCH --account ELIX5_fratelli
#SBATCH --error wgcna_ERR
#SBATCH --output wgcna_OUT

# Modules and environment activation (bash)
module load profile/bioinf
module load r/4.1.0--gcc--10.2.0-python--3.8.6

# Funzione principali
Rscript WGCNA_mRNA/1_create_object.R

exit