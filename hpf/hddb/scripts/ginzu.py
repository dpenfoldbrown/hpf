#!/usr/bin/env python

import os
import sys
import getopt
import subprocess
from hpf.runtime import runtime, Flag, IntegerOption, FileOption, Option
from hpf.utilities import consume
from hpf.processing import processor, SYNCHRONOUS
from hpf.hddb import proteins, connect, tunnel
tunnel(10)
from hpf.hddb.db import Protein, Experiment, Session

SEQUENCES = "sequences"
ORGANISM_TYPE = "organism_type"
GINZU_VERSION = "ginzu_version"

def ginzu(sequence_key,sleep=10):
    # Randomize the calls to tunnel so we don't flood with multiple connections
    cmd = "ddb.pl -mode sequence -submode process -sequence_key {0} -ginzu_version {1}".format(sequence_key, runtime().opt(GINZU_VERSION))
    if runtime().opt(ORGANISM_TYPE):
        cmd += " -organismtype %s" % runtime().opt(ORGANISM_TYPE)
    runtime().debug(cmd)
    retcode = subprocess.check_call(cmd, cwd=os.getcwd(), shell=True)
    return sequence_key
    
def tasks(*experiment_key):
    session = Session()
    proteins = []
    if len(experiment_key)>0:
        proteins += session.query(Protein).join(Protein.experiment).filter(Experiment.id.in_(map(int,experiment_key))).distinct().all()
        runtime().debug("Found %i proteins in experiments" % len(proteins),experiment_key)
        #DEBUG#
        print "Found {0} proteins in experiments {1}".format(len(proteins),experiment_key) 
    if runtime().opt(SEQUENCES) != None:
        with open(runtime().opt(SEQUENCES)) as handle:
            keys = [int(i.strip()) for i in handle.readlines()]
            runtime().debug("Found %i keys in file" % len(keys))
            proteins += session.query(Protein).filter(Protein.sequence_key.in_(keys)).distinct().all()
    return list(set([protein.sequence_key for protein in proteins]))
#    
#    with connect() as cursor:
#        query = """select distinct sequence_key from hpf.protein where experiment_key in (%s)""" % ",".join(experiment_key)
#        runtime().debug(query)
#        cursor.execute(query)
#        return [r[0] for r in cursor.fetchall()]
#        #return [protein.id for protein in proteins(cursor, experiment=experiment_key, filter_experiments=False)] 

def main(*experiment_key):
    # Check for ginzu_version. If no, die.
    if not runtime().opt(GINZU_VERSION):
        raise "Fatal exception: must provide valid -g GINZU_VERSION argument"
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS))
    runtime().debug("Using processor",pool)
    pool.make_tasks(tasks,*experiment_key)
    consume(pool.run(ginzu))
    
def _do(argv):
    r = runtime()
    r.description("""
    ginzu.py [-options] *experiment_key
    Run ginzu for any number of experiment keys.
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing."))
    r.add_option(FileOption(SEQUENCES, "f", description="Integer on each line, sequence keys to process.",default=None))
    r.add_option(Option(ORGANISM_TYPE, "o", description="Type of organism for SignalP runs (gram+, gram-, euk)", default=""))
    r.add_option(Option(GINZU_VERSION, "g", description="REQUIRED. Ginzu version to run ginzu under. Must be in DB, or ginzu fails",default=None))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
