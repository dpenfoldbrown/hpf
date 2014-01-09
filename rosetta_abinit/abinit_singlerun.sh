#!/bin/bash

# A simple bash driver to run Rosetta 3.3 Ab init structure prediction
# Creates a flags file then runs, as it is nice to have a record of
# the params you used to generated struct preds (ie, that flags file)
# NOTE as always that PSIPRED file must be psipred ss2 format

ABINIT_EXEC="/share/apps/rosetta/RosettaReleases/rosetta-3.3/rosetta_source/bin/AbinitioRelax.default.linuxgccrelease"
ROSETTA_DB="/share/apps/rosetta/RosettaReleases/rosetta-3.3/rosetta_database"

WORK_DIR="/home/dpb3/rosetta/work/qv096"
FLAGS_FILE="$WORK_DIR/abinit_flags"

FASTA_IN="$WORK_DIR/qv096.fasta"
FRAG3_IN="$WORK_DIR/aaqv09603_05.200_v1_3"
FRAG9_IN="$WORK_DIR/aaqv09609_05.200_v1_3"
PSIPRED_IN="$WORK_DIR/qv096.ss2"

# Assuming using silent file out format (for PDB, change structure)
OUTFILE="$WORK_DIR/qv096_silent.out"

# Number of decoy structures (models) to produce
MODELS="3"


## Write flags file

echo "-in:file:fasta $FASTA_IN"       > $FLAGS_FILE
echo "-in:file:frag3 $FRAG3_IN"      >> $FLAGS_FILE
echo "-in:file:frag9 $FRAG9_IN"      >> $FLAGS_FILE
echo "-database $ROSETTA_DB"         >> $FLAGS_FILE
echo ""                              >> $FLAGS_FILE
echo "-abinitio:relax"               >> $FLAGS_FILE
echo "-relax:fast"                   >> $FLAGS_FILE
echo "-abinitio::increase_cycles 10" >> $FLAGS_FILE
echo "-abinitio::rg_reweight 0.5"    >> $FLAGS_FILE
echo "-abinitio::rsd_wt_helix 0.5"   >> $FLAGS_FILE
echo "-abinitio::rsd_wt_loop 0.5"    >> $FLAGS_FILE
echo ""                              >> $FLAGS_FILE
echo "-use_filters true"             >> $FLAGS_FILE
echo "-psipred_ss2 $PSIPRED_IN"      >> $FLAGS_FILE
echo "-kill_hairpins $PSIPRED_IN"    >> $FLAGS_FILE
echo ""                              >> $FLAGS_FILE
echo "-out:file:silent $OUTFILE"     >> $FLAGS_FILE
echo "-nstruct $MODELS"              >> $FLAGS_FILE


## Run Rosetta ab init struct pred

cd $WORK_DIR
echo "Start: " `date`  > abinit_time
$ABINIT_EXEC @$FLAGS_FILE
EXIT=$?
echo "Finis: " `date` >> abinit_time
exit $EXIT

