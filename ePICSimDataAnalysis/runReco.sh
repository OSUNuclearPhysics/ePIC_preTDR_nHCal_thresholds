#!/bin/bash

#/opt/benchmarks/physics_benchmarks/bin/get_calibrations

#export DETECTOR_PATH=/gpfs/mnt/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev/epic/install/share/epic


export LOCAL_PREFIX=/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev

source ${LOCAL_PREFIX}/ip6/install/setup.sh
source ${LOCAL_PREFIX}/epic/install/setup.sh

#export DETECTOR_CONFIG=epic_full

export JUGGLER_SIM_FILE=data/nhcal_sim.edm4hep.root
export JUGGLER_REC_FILE=output/nhcal_sim_neutron.reco.edm4hep.root
export JUGGLER_N_EVENTS=200


#gaudirun.py /opt/benchmarks/physics_benchmarks/options/reconstruction.py

gaudirun.py reconstruction.py