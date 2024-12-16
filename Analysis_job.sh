#!/bin/bash
#SBATCH --job-name=analysis_job
#SBATCH --account=PAS2524
#SBATCH --time=00:10:00

/users/PAS2524/corey90/ePIC_preTDR_nHCal_thresholds/container/eic-shell <<EOF
/users/PAS2524/corey90/ePIC_preTDR_nHCal_thresholds/ePICSimDataAnalysis/readTreeSimMain ${INPUT_LIST} ${OUTPUT_FILE} 80000
EOF
