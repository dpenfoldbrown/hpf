#!/bin/bash

# A simple driver for Rosetta 3.3 fragment picker via make_local_fragments.pl

FRAG_SCRIPT="/home/dpb3/rosetta/rosetta_fragments/make_fragments.pl"
FLAGS="-verbose -nosam -noporter -nocleanup -nohoms"

OUT_DIR="/home/dpb3/rosetta/work/rc097"
PSIPRED_IN="$OUT_DIR/rc097.psipred.2" 
FASTA_IN="$OUT_DIR/rc097.fasta"

$FRAG_SCRIPT $FLAGS -psipredfile $PSIPRED_IN -rundir $OUT_DIR $FASTA_IN

