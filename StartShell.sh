#!/bin/bash

OUT_FILE=$1
ENERGY=$2
N_EVT=$3
LIST_FILE=$4

$PROJECT_DIR/container/eic-shell <<EOF
source $PROJECT_DIR/epic/install/bin/thisepic.sh
npsim --compactFile=$PROJECT_DIR/epic/install/share/epic/epic_backward_hcal_only.xml -N=$N_EVT --enableGun --gun.particle=\"neutron\" --gun.energy $ENERGY*GeV --gun.thetaMin 130*deg --gun.thetaMax 177*deg --gun.distribution uniform --outputFile $OUT_FILE
echo "$OUT_FILE" >> $LIST_FILE
EOF
