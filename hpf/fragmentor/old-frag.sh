#!/bin/bash
echo "HOSTNAME = $HOSTNAME"
echo "SHELL = $SHELL"
echo "PATHS = $PATH"
echo "CWD = $PWD"
echo "PYTHONPATH = $PYTHONPATH"
echo "+++++++++++++++++++++++++++++++++++++++++"

# Script arguments: frag_array.py <prediction code> <fasta file input> <output location>
# prediction code is a two letter code, such as 'mm'
# Note: the %s in the fasta file input will be replaced by the prediction code and task ID, e.g: mm078

# To run, set OUTFILECODE env to a valid code (eg, oi), and run/qsub the script.

OUTFILECODE="oi"
TEST_STR="_test"
INPUT="/scratch/fragmentor/input/$OUTFILECODE$TEST_STR/%s.fasta"
OUTPUT_DIR="/scratch/fragmentor/results/old"

echo "Outfile code     : $OUTFILECODE"
echo "Input string     : $INPUT"
echo "Output directory : $OUTPUT_DIR"

python ./src/old/frag_array.py $OUTFILECODE $INPUT $OUTPUT_DIR
exit
