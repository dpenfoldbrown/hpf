#!/usr/bin/env python

# A driver for parallelizing interpro runs over the cluster via the HPF task pool
# structure. Creates an Interpro object from hpf.interpro.interpro.py

from hpf.utilities import consume
from hpf.processing import processor
from hpf.hddb.db import Session, Protein
from hpf.interpro.interpro import Interpro

INTERPRO_DIR = "/scratch/interpro"

def interpro_driver(sequence_id):
    print "Interpro driver::interpro_driver:: Creating and running Interpro object on sequence {0}".format(sequence_id)
    ipr = Interpro(sequence_id=sequence_id, results_dir=INTERPRO_DIR+"/interpro_results")
    try:
        ipr.run()
    except:
        print "Interpro run on sequence {0} failed".format(sequence_id)
        raise

def pfam_driver(sequence_id):
    print "Interpro driver::pfam_driver:: Creating and running Interpro-pfam object on sequence {0}".format(sequence_id)
    pfam = Interpro(sequence_id=sequence_id, results_dir=INTERPRO_DIR+"/pfam_results", appl_str="hmmpfam")
    try:
        pfam.run()
    except:
        print "Interpro-pfam run on sequence {0} failed".format(sequence_id)
        raise


def tasks(experiment_id):
    print "Interpro driver::tasks:: Getting tasks (sequences) for experiment {0}".format(experiment_id)
    session = Session()
    proteins = session.query(Protein.sequence_key).filter_by(experiment_key=experiment_id)
    experiment_seqs = list()
    for protein in proteins:
        experiment_seqs.append(protein[0])
    if proteins.count() != len(experiment_seqs):
        raise Exception("Number of sequences extracted from experiment {0} does not match number of proteins".format(experiment_id))
    unique_seqs = list(set(experiment_seqs))
    print "{0} tasks (unique sequences) from {1} sequences retrieved".format(len(unique_seqs), len(experiment_seqs))
    return unique_seqs


def main():
    # TODO: make a pfam flag commandline argument
    # TODO: make output dir a commandline argument
    # TODO: make experiment a commandline argument
    experiment_id = 1171

    print "Interpro driver:: Creating processor pool, no SYNCHRONOUS option"
    pool = processor()
    pool.make_tasks(tasks, experiment_id)
    print "Running tasks..."
    #consume(pool.run(interpro_driver))
    consume(pool.run(pfam_driver))


if __name__ == "__main__":
    main()

