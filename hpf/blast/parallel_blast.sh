#!/bin/bash
# Driver script for executing parallel_blast.py

echo "Hostname = $HOSTNAME"
echo "Shell = $SHELL"
echo "Python path: $PYTHONPATH"

date
echo "--- PARALLEL BLAST ---"
python /home/dpb3/hpf/blast/parallel_blast.py
echo "--- Complete       ---"
date
exit $?

