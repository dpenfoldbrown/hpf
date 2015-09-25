#!/bin/bash
echo "HOST: " `hostname`
echo "SHELL: $SHELL"
echo "Python path: $PYTHONPATH"
echo "BLAST: " `which blastpgp`
echo "PSIPRED: " `which psipred`
echo "mcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcmmcm"


DRIVER="/home/dpb3/hpf/hpf/mcm/mcmrun_driver.py"

CODE="ox"
VERSION=1
HOST="HOSTNAME"
WORK_DIR="/scratch/dpb3/mcm/work"
RESULTS_DIR="/scratch/dpb3/mcm/results/ox_hpf2_results"


echo "Code       : $CODE"
echo "Version    : $VERSION"
echo "Work dir   : $WORK_DIR"
echo "Results dir: $RESULTS_DIR"

date

python $DRIVER --code $CODE --version $VERSION --host $HOST --work_dir $WORK_DIR --results_dir $RESULTS_DIR

#
# To re-run error codes (add flag)
# python $DRIVER --code $CODE --version $VERSION --work_dir $WORK_DIR --results_dir $RESULTS_DIR --redo_error
#

EXIT=$?
date
exit $EXIT
