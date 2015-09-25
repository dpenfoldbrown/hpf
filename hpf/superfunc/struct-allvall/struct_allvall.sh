#!/bin/bash
# Driver script for executing homolog_struct_allvall.py

echo "Hostname = $HOSTNAME"
echo "Shell = $SHELL"
echo "Python path: $PYTHONPATH"

date
echo "--- HOMOLOG STRUCTUR ALL v ALL ---"
python /home/dpb3/hpf/superfunc/cluster_struct_allvall.py
echo "--- Complete       ---"
date
exit $?
