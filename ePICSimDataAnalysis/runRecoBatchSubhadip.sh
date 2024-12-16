#!/bin/sh

source /opt/detector/epic-nightly/setup.sh 
source /gpfs/mnt/gpfs02/eic/palsp/eicsoft/epic/install/setup.sh 

gaudirun.py /opt/benchmarks/physics_benchmarks/options/reconstruction.py 
 
#EOF