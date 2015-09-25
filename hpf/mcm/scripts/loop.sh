#!/bin/bash
for ((i=$2;i<=$3;i++)); do
qsub -vCODE=$1,TASK=$i -N mcm_$1 array.sh
done
