#!/usr/bin/env python

import sys, os
import getopt
from hpf.runtime import runtime, Flag
from hpf.utilities import consume
from hpf.hddb.db import *
from Bio import SeqIO
from hashlib import sha1

session = None

from hashlib import sha1
def digest(record):
    return sha1(str(record.seq)).hexdigest()

def map(function, sequence):
    i=0
    for item in sequence:
        yield function(item)
        i+=1
        if i%100 == 0:
            print "finished %i" % i

def rename(record):
    #print record.id
    try:
        s = session.query(Sequence).join(Sequence.ac).filter(SequenceAc.ac==record.id).one()
    except:
        print "Can't find ",record.id
        try:
            s = session.query(Sequence).filter(Sequence.sha1==digest(record)).one()
        except:
            print record.id
            print record.seq
    record.id=str(s.id)
    return record

        
def main(*files):
    global session
    #runtime().set_debug(1)
    session = Session()
    for file in files:
        runtime().debug("Mapping file",file)
        with open(file) as handle:
            records = list(SeqIO.parse(handle,"fasta"))
        runtime().debug("Found %i records" % len(records))
        with open(file+".hpf","w") as handle:
            SeqIO.write(map(rename,records), handle, "fasta")
            
def _do(argv):
    r = runtime()
    r.description("""
    map_fasta.py [-options] fasta
    Rename all of the records with sequence id's from the database.
    """)
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
