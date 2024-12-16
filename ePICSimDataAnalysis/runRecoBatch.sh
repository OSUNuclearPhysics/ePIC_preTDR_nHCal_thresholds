#!/bin/bash

#curl -L get.athena-eic.org | bash


EICSHELL=/gpfs/mnt/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev/eic-shell
#EICSHELL=./eic-shell

#source /opt/detector/epic-nightly/setup.sh

export LOCAL_PREFIX=/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev

source ${LOCAL_PREFIX}/ip6/install/setup.sh
source ${LOCAL_PREFIX}/epic/install/setup.sh

#export DETECTOR=epic
#export DETECTOR_PATH=/gpfs/mnt/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev/epic/install/share/epic
#export DETECTOR_PATH=/gpfs/mnt/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev/local/share/epic
#export ATHENA_PREFIX=/gpfs/mnt/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev/local
#export DETECTOR_PREFIX=/gpfs/mnt/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev//epic/install
#export DETECTOR_PATH=local/share/epic
#export ATHENA_PREFIX=local
#export DETECTOR_VERSION=master
#export BEAMLINE_CONFIG=ip6
#export BEAMLINE_CONFIG_VERSION=master

## note: we will phase out the JUGGLER_* flavor of variables in the future
export JUGGLER_DETECTOR=$DETECTOR
export JUGGLER_DETECTOR_VERSION=$DETECTOR_VERSION
export JUGGLER_DETECTOR_PATH=$DETECTOR_PATH
export JUGGLER_INSTALL_PREFIX=${ATHENA_PREFIX}
export JUGGLER_BEAMLINE_CONFIG=$BEAMLINE_CONFIG
export JUGGLER_BEAMLINE_CONFIG_VERSION=$BEAMLINE_CONFIG_VERSION

## Export detector libraries
#export LD_LIBRARY_PATH=$DETECTOR_PREFIX/lib:$LD_LIBRARY_PATH

## modify PS1 for this detector version
export PS1="${PS1:-}"
export PS1="nightly${PS1_SIGIL}>${PS1#*>}"
unset branch

CONDOR_DIR=condorReco
OUT_DIR=output

mkdir ${CONDOR_DIR}
mkdir "${CONDOR_DIR}/${OUT_DIR}"

cp -r local/share/ip6 .

#cat << EOF | ${EICSHELL}

#${EICSHELL} <<EOF

#/opt/benchmarks/physics_benchmarks/bin/get_calibrations

export JUGGLER_SIM_FILE=${1}
export JUGGLER_REC_FILE=${2}
export JUGGLER_N_EVENTS=${3}

#echo "gaudirun.py reconstruction.py" | ${EICSHELL} | tee "/gpfs/mnt/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev/condorReco/${5}_${4}.log"
echo "gaudirun.py reconstruction.py" | ${EICSHELL}

#EOF