#!/bin/bash

# A run script for Rosetta 3.3 ab initio structure prediction
# designed to be executed on cluster, where many runs on the same
# sequence are happening simultaneously

# Requires a pre-existing flags file (with -out UNSPECIFIED)

ABINIT_EXEC="/share/apps/rosetta/RosettaReleases/rosetta-3.3/rosetta_source/bin/AbinitioRelax.default.linuxgccrelease"

WORK_DIR="/home/dpb3/rosetta/work/rc097"
OUT_DIR="$WORK_DIR/abinit_results"

CODE="rc097"
FLAGS_FILE="$WORK_DIR/abinit_flags"
OUT_FILE="$OUT_DIR/$CODE.silent.$PBS_JOBID"


# Check immediately for PBD JOB ID. Can't work without it
if [ -z "$PBS_JOBID" ]; then
    echo "PBS_JOBID not present, can not uniquely identify run"
    exit 1
fi

# Check for flags file
if [ ! -f "$FLAGS_FILE" ]; then
    echo "Flags file $FLAGS_FILE is not valid (DNE)"
    exit 2
fi

# Check/make silent file output directory in work dir
if [ ! -d "$OUT_DIR" ]; then
    mkdir $OUT_DIR
fi


# Messages
echo "Running ab init struct pred on code $CODE"
echo "Run information in flags at $FLAGS_FILE"
echo "Job $PBS_JOBID"
echo "Ab init outfile: $OUT_FILE"


# Run Rosetta ab init struct pred
cd $WORK_DIR
echo "Start: " `date`
$ABINIT_EXEC @$FLAGS_FILE -out:file:silent $OUT_FILE
EXIT=$?
echo "Finis: " `date`
exit $EXIT

