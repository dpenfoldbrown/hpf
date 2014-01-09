#!/bin/bash
out=`ps -ef | grep "ssh -i /home/kd35/.ssh/id_rsa_usq -f -N -L 13307:mysql:3306 kdrew@db.proteomics.washington.edu" | wc -l`
#echo $out
if [ $out -gt 1 ]; then 
    pkill -o -f "ssh -i /home/kd35/.ssh/id_rsa_usq -f -N -L 13307:mysql:3306 kdrew@db.proteomics.washington.edu"
fi
exit $?
