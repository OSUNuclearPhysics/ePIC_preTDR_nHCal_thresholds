#!/bin/bash
#SBATCH --job-name=npsim_job
#SBATCH --account=PAS2524
#SBATCH --time=02:00:00
#SBATCH --array=1-50

PROJECT_DIR="/users/PAS2524/corey90/ePIC_preTDR_nHCal_thresholds"
export PROJECT_DIR

#$PROJECT_DIR/StartShell.sh ${OUTPUT_DIR}/neutron_${ENERGY}GeV_${SLURM_ARRAY_TASK_ID}.root $ENERGY $N_EVT $LIST_FILE
$PROJECT_DIR/StartShell.sh "${OUTPUT_DIR}/neutron_${ENERGY}GeV_${SLURM_ARRAY_TASK_ID}.root" "$ENERGY" "$N_EVT" "$LIST_FILE"
