#!/bin/bash

PROJECT_DIR="/users/PAS2524/corey90/ePIC_preTDR_nHCal_thresholds"

OUTPUT_DIR="${PROJECT_DIR}/NPSim_outputs"
mkdir -p "$OUTPUT_DIR"
LIST_DIR="${PROJECT_DIR}/ListFiles"
mkdir -p "$LIST_DIR"

N_EVT=2000

# Define the energies you want to simulate
ENERGIES=(0.1 0.2 0.3 0.4 0.5 0.65 0.8 1 1.5 2 3 4 5)

for ENERGY_VAL in "${ENERGIES[@]}"; do
    # Create the list file
    LIST_FILE="${LIST_DIR}/bhcal_${ENERGY_VAL}gev.list"
    touch "$LIST_FILE"

    # Submit the job
    sbatch --export=OUTPUT_DIR="$OUTPUT_DIR",ENERGY="$ENERGY_VAL",N_EVT="$N_EVT",LIST_FILE="$LIST_FILE" "$PROJECT_DIR/SubmitFile.sh"
    #$PROJECT_DIR/SubmitFile.sh "$OUTPUT_DIR" "$ENERGY_VAL" "$N_EVT" "$LIST_FILE"
done
