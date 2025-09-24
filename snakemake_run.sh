#!/bin/bash

## Conda environment
conda activate BeeGees_env
 
# Setup logging
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
VERSION="v1.0.0"
RUN_ID="BeeGees Snakemake workflow"				
LOG_FILE="snakemake_${TIMESTAMP}.log"

# Function for timestamped logging
log_with_timestamp() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Start logging
log_with_timestamp "=== Snakemake Workflow Start ==="
log_with_timestamp "Version: $VERSION"
log_with_timestamp "Run ID: $RUN_ID"
log_with_timestamp "Log file: $LOG_FILE"
log_with_timestamp "Screen session: $STY"

# Show system info
log_with_timestamp "Running on: $(hostname)"
log_with_timestamp "Working directory: $(pwd)"
log_with_timestamp "Conda environment: $CONDA_DEFAULT_ENV"

##Snakemake
# Unlock directory
log_with_timestamp "Unlocking Snakemake directory..."
snakemake --profile ./profiles/slurm/ \
    --snakefile ./workflow/Snakefile \
	--configfile ./config/config.yaml \
	--unlock

# Run snakemake workflow with profile
log_with_timestamp "Starting workflow execution..."
snakemake --profile ./profiles/slurm/ \
    --snakefile ./workflow/Snakefile \
    --configfile ./config/config.yaml \
    --rerun-incomplete \
    2>&1 | tee -a "$LOG_FILE"

log_with_timestamp "Log file saved as: $LOG_FILE"
