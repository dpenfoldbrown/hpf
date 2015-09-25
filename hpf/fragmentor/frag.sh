#!/bin/bash
echo "HOSTNAME = " `hostname`
echo "SHELL = $SHELL"
echo "PYTHONPATH = $PYTHONPATH"
echo "fragmentorfragmentorfragmentorfragmentorfragmentorfragmentorfragmentorfragmentor"

# Script arguments: frag_driver.py <prediction code> <fasta file input> <output location>
# prediction code is a two letter code, such as 'mm'

OUTFILECODE="qr"
INPUT="/scratch/dpb3/fragmentor/input/$OUTFILECODE/*.fasta"
OUTPUT_DIR="/scratch/dpb3/fragmentor/results/"

echo "Outfile code     : $OUTFILECODE"
echo "Input string     : $INPUT"
echo "Output directory : $OUTPUT_DIR"

python ./frag_driver.py $OUTFILECODE "$INPUT" $OUTPUT_DIR
exit $?
