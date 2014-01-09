#!/usr/bin/python

# Fetches sequence ID, Accession, and Gene Name (if found in sequenceAc.description)
# from the database, given hpf experiment ID

import re
import sys
import argparse

def main(experiment, outfile):
    from hpf.hddb.db import Session, Protein
    session = Session()
    
    gn_pattern = r"GN=(?P<gene_name>\S+)"
    
    handle = open(outfile, 'w')
    proteins = session.query(Protein).filter_by(experiment_key=experiment)
    
    for p in proteins:
        gn_found = re.search(gn_pattern, p.ac.description)
        gene_name = gn_found.group('gene_name') if gn_found else "None"
        #if p.ac.ac2:
        #    handle.write("{0}\t{1}\t{2}\n".format(p.sequence_key, p.ac.ac2, gene_name))
        #else:
        handle.write("{0}\t{1}\t{2}\n".format(p.sequence_key, p.ac.ac, gene_name))
        sys.stdout.write("protein: {0}    {1}               \r".format(p.sequence_key, gene_name))
    
    handle.close()
    print "Map complete"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A tool to map HPF sequence key proteins for an experiment to AC and Gene name")
    parser.add_argument("--experiment", action="store", dest="experiment_key", required=True, type=int,
                help="Experiment key of proteins to create map for")
    parser.add_argument("--outfile", action="store", dest="outfile", default=None,
                help="File to write map to")
    args = parser.parse_args()

    if not args.outfile:
        args.outfile = "{0}_hpf_ac_gn.txt".format(args.experiment_key)
    main(args.experiment_key, args.outfile)
