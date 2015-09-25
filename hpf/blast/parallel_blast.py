#!/usr/bin/env python

# A cluster parallelization script for blasting on the cluster.
# Runs on a filestring-specified collection of blast query-sequence fasta files
# (eg: 'somedir/*.fasta' will run an individual blast on the cluster for each
# (file that *.fasta expands to)
# Not as clean or good as it should be, time constrained.
#
# See the split-fasta script to split a large fasta file into smaller pieces
# dpb 10/25/2011

import sys
from glob import glob
from os.path import expanduser, isfile
from hpf.utilities import consume, system
from hpf.processing import processor

# Blast paramaters
blast_exe = 'blastpgp'
source_db = '/scratch/jamboree/databases/fourorgs_doms.fasta'
e_val = '1e-3'
results = 250
alignments = 250
processors = 1
output_mode= 7


def run_blast(query_file):
    if not isfile(query_file):
        raise Exception("Given blast query file {0} not a valid file".format(query_file))
    blast_outfile = query_file + '.blast.xml'
    blast_cmd = "{0} -i {1} -d {2} -e {3} -v {4} -b {5} -a {6} -m {7} -o {8}".format( \
                blast_exe, query_file, source_db, e_val, results, alignments, processors, output_mode, blast_outfile)
    print "Blast command: {0}".format(blast_cmd)
    try:
        system(blast_cmd)
    except Exception as e:
        print e
        print "Blasting query file {0} against DB {1} failed".format(query_file, source_db)
        raise


def tasks(queryfstr):
    query_files = glob(expanduser(queryfstr))
    if query_files == []:
        raise Exception("Input string {0} expands to no valid files".format(queryfstr))
    print "tasks:: input blast query files: ", query_files
    print "tasks:: {0} tasks set".format(len(query_files))
    return query_files


def main():

    queryfile_filestring = '/scratch/jamboree/fasta_split/hpd.fasta.*'

    pool = processor()
    pool.make_tasks(tasks, queryfile_filestring)
    consume(pool.run(run_blast))


def test():
    query_files = ['/home/dpb3/superfunc/data/split_human1176_500/human_1176_c90s80_RNA.fasta.156',
        '/home/dpb3/superfunc/data/split_human1176_500/human_1176_c90s80_RNA.fasta.166',
        '/home/dpb3/superfunc/data/split_human1176_500/human_1176_c90s80_RNA.fasta.343'
        ]
    
    for file in query_files:
        run_blast(file)
    print "Run blast parallel test complete"


if __name__ == "__main__":
    #test()
    main()

