#!/bin/bash

# A quick GINZU single-run script
# ddb.pl must be on PATH to run

SEQUENCE_KEY=17397097
GINZU_VERSION=4

date
ddb.pl -mode sequence -submode process -sequence_key $SEQUENCE_KEY -ginzu_version $GINZU_VERSION
exit $?
