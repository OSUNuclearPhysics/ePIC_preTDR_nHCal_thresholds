#!/bin/bash

PROJECT_DIR="/users/PAS2524/corey90/ePIC_preTDR_nHCal_thresholds"
LIST_DIR="${PROJECT_DIR}/ListFiles"
OUTPUT_DIR="${PROJECT_DIR}"

# Define the energies you want to analyze
ENERGIES=(0.1 0.2 0.3 0.4 0.5 0.65 0.8 1 1.5 2 3 4 5)

for ENERGY_VAL in "${ENERGIES[@]}"; do
    # Define input and output files
    INPUT_LIST="${LIST_DIR}/bhcal_${ENERGY_VAL}gev.list"
    OUTPUT_FILE="${OUTPUT_DIR}/backhcal_${ENERGY_VAL}gev_batch1.root"

    # Submit the job
    sbatch --export=INPUT_LIST="$INPUT_LIST",OUTPUT_FILE="$OUTPUT_FILE" "$PROJECT_DIR/Analysis_job.sh"
done
