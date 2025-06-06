#!/bin/bash
#SBATCH --job-name=BGEE-controller
#SBATCH --partition=long


## Conda environment
source ~/conda/etc/profile.d/conda.sh

conda activate bgee_env


##Snakemake
# Unlock working directory
snakemake --profile ./profiles/slurm/ --snakefile ./workflow/Snakefile-196beta3 --configfile ./config/config.yaml --unlock


# Run snakemake workflow with profile
snakemake --profile ./profiles/slurm/ \
    --snakefile ./workflow/Snakefile-196beta3 \
    --configfile ./config/config.yaml
