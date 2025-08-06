#!/bin/bash

## Activate conda environment
source path/to/conda/etc/profile.d/conda.sh
conda activate bgee_env


# Set run ID
RUN_ID="BGEE_RUN_1"	

# Set Snakemake profile (slurm or local)
PROFILE=".profiles/slurm/"


# Setup logging
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="snakemake_${TIMESTAMP}.log"

log_with_timestamp() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

log_with_timestamp "=== Snakemake Workflow Start ==="
log_with_timestamp "Run ID: $RUN_ID"
log_with_timestamp "Log file: $LOG_FILE"
log_with_timestamp "Running on: $(hostname)"
log_with_timestamp "Working directory: $(pwd)"
log_with_timestamp "Conda environment: $CONDA_DEFAULT_ENV"




## Snakemake
# Unlock directory
log_with_timestamp "Unlocking Snakemake directory..."

snakemake --profile $PROFILE --snakefile ./workflow/Snakefile --configfile ./config/config.yaml --unlock


# Run snakemake workflow with profile
log_with_timestamp "Starting workflow execution..."

snakemake --profile $PROFILE \
    --snakefile ./workflow/Snakefile \
    --configfile ./config/config.yaml \
    --rerun-incomplete \
    --cores all \
    2>&1 | tee -a "$LOG_FILE"

log_with_timestamp "Log file saved as: $LOG_FILE"
