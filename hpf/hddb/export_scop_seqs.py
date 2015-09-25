#!/usr/bin/env python

# A script to export SCOP sequences from the HPF database to file in fasta 
# format. hpf.astral table (db.Astral) contains sccs and sequence key for all
# SCOP SFs

import argparse

def main():
    from hpf.hddb.db import Session, Sequence, Astral
    session = Session()
    astral_records = session.query(Astral)
    with open(args.outfile, 'w') as handle:
        for astral in astral_records:
            print "{0}\r".format(astral.sccs),
            handle.write(">{0}|{1}{2}|{3}\n".format(astral.sccs, astral.pdbid, astral.chain, astral.sequence_key))
            handle.write("{0}\n".format(astral.sequence.sequence))
    print "Exporting SCOP sequences to file {0} complete".format(args.outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to export SCOP sequences from HPF DB")
    parser.add_argument("-o", "--outfile", action="store", dest='outfile', required=True,
            help='File to output SCOP fasta records to')
    args = parser.parse_args()
    main()
