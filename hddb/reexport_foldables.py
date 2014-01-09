#!/usr/bin/env python

# Script to re-export foldables already run through the ddb.pl -mode export -submode foldables
# program. Works by two-letter code (eg, 'qa') on the command line.
# Exports the foldable sequence from the filesystemOutfile table (ORM: FilesytemOutfile)
#
# NOTE: FOLDABLE SEQUENCE IS NOT THE DOMAIN SEQUENCE. A foldable domain sequence is pared 
# down (remove signal sequence, transmembrane, etc) to create the foldable sequence
#
# dpb, 8/08/2012

import os
import re
import argparse

def _init_db():
    """Initialize a database session and return session object"""
    from hpf.hddb.db import Session
    print "Initializing database session"
    return Session()

def write_fasta(handle, seqid, length, sequence):
    """Writes a given sequence and corr. information to the given handle in a specific format"""
    handle.write(">sequence.id.{0} (database: hpf; length: {1};)\n{2}\n".format(seqid, length, sequence))


def export_single(prediction_code):
    """Exports a fasta file with the sequence corresponding to given code"""
    from hpf.hddb.db import FilesystemOutfile
    session = _init_db()
    foldable_record = session.query(FilesystemOutfile).filter_by(prediction_code=prediction_code).first()
    if not foldable_record:
        raise Exception("Foldable record with code {0} not found in the database".format(prediction_code))
    with open(prediction_code+".fasta", 'w') as handle:
        write_fasta(handle, foldable_record.sequence_key, len(foldable_record.sequence.sequence), foldable_record.sequence.sequence)
    print "Fasta {0} written for code {1}".format(prediction_code+".fasta", prediction_code)


def export_all(prediction_code):
    """Exports fasta files for all code of the letter code. EG: aa => aa000, aa001... aa999"""
    from hpf.hddb.db import FilesystemOutfile
    session = _init_db()
    foldables = session.query(FilesystemOutfile).filter(FilesystemOutfile.prediction_code.like(prediction_code+'%')).all()
    print "{0} foldable records found for code {1}".format(len(foldables), prediction_code)
    written = 0
    for fr in foldables:
        with open(fr.prediction_code+".fasta", 'w') as handle:
            write_fasta(handle, fr.sequence_key, len(fr.sequence.sequence), fr.sequence.sequence)
        print "Fasta {0} written for code {1}".format(fr.prediction_code+".fasta", fr.prediction_code)
        written += 1
    print "{0} fasta files written for letter code {1}".format(written, prediction_code)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A tool to re-export already processed foldable sequences")
    parser.add_argument("-m", "--mode", action="store", dest="mode", default="all", choices=["single", "all"],
            help="Single exports one fasta only (code must be whole, eg aa999). All exports all fasta for a two-letter code")
    parser.add_argument("-c", "--prediction_code", action="store", dest="prediction_code", required=True,
            help="The prediction code to export fasta files for. Single mode: 'aa999'. All mode: 'aa'")
    parser.add_argument("-d", "--directory", action="store", dest="directory", default=os.getcwd(),
            help="Directory to export fasta files to")
    args = parser.parse_args()
    
    
    if not os.path.isdir(args.directory):
        raise EnvironmentError("Directory {0} is not a valid directory".format(args.directory))
    os.chdir(args.directory)

    if args.mode == 'single':
        if not re.match(r"[a-z]{2}[0-9]{3}", args.prediction_code):
            raise Exception("Code {0} does not match prediction code format, EG: aa999".format(args.prediction_code))
        export_single(args.prediction_code)
    
    elif args.mode == 'all':
        if not re.match(r"[a-z]{2}", args.prediction_code):
            raise Exception("Code {0} does not match two-letter code format, EG: aa".format(args.prediction_code))
        export_all(args.prediction_code)

