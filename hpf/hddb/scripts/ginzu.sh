#!/bin/bash
echo "Hostname = $HOSTNAME"
echo "Shell = $SHELL"
echo "Python path: $PYTHONPATH"
echo "ginzuginzuginzuginzuginzuginzuginzuginzuginzuginzuginzuginzuginzuginzu"

# Script to drive the execution of ginzu.py (qsub often can not take script arguments).
# script arguments: ginzu.py <experiment_code>

GINZU="/home/dpb3/ginzu/hddb/scripts/ginzu.py"
GINZU_VERSION="4"
EXPERIMENT="1179"

echo "Execution directory: "
pwd 
echo "Experiment in DB to ginzu: $EXPERIMENT"
echo "Ginzu version given: $GINZU_VERSION"

date
python $GINZU -g $GINZU_VERSION $EXPERIMENT
date
exit $?

