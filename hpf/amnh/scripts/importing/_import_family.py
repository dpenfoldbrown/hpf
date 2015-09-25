#!/usr/bin/env python

from hpf.hddb.db import *
import sys, os
import getopt
from hpf.runtime import runtime, Flag, debug
from hpf.utilities import consume
from hpf.processing import processor, SYNCHRONOUS
import time
from Bio import SeqIO
import hashlib
import string
from sqlalchemy.sql import *

E_ID = 1162 # Hard coded to match temporary AMNH experiment
Session = None

def _import_family(file):
    thread_setup()
    return import_family(file)

def thread_setup():
    global Session
    if Session==None:
        clear()
        e,b,s = setup()
        debug("Entering subprocess setup",s,Session)
        Session = s
    

def import_family(file):
    name = ".".join(os.path.basename(file).split(".")[:-1]).replace("_"," ")
    debug(name, file)
    session = Session()
    family = session.query(Family).filter(Family.name==name).first()
    if family==None:
        family = Family()
        family.name = name
        family.experiment_key = E_ID
        session.add(family)
        session.flush()
    
    with open(file) as handle:
        for record in SeqIO.parse(handle,"fasta"):
            seq = str(record.seq).replace("*","")
            if not all([c in string.ascii_uppercase for c in seq]):
                runtime().println("Malformed Sequence", Family, seq)
            sha1 = hashlib.sha1(seq).hexdigest()
            sequence = session.query(Sequence).filter(Sequence.sha1==sha1).first()
            if sequence==None:
                sequence = Sequence()
                sequence.sha1 = sha1
                sequence.sequence=seq
                id = session.add(sequence)
                session.flush()
                debug("Added",sequence)
            sequenceAc = session.query(SequenceAc).filter(and_(SequenceAc.sequence_key==sequence.id, SequenceAc.ac==record.id)).first()
            if sequenceAc==None:
                sequenceAc = SequenceAc()
                sequenceAc.sequence_key=sequence.id
                sequenceAc.gi = None
                sequenceAc.db = "amnh"
                sequenceAc.ac = record.id
                sequenceAc.description = "amnh"
                sequenceAc.taxonomy_id = 0
                session.add(sequenceAc)
                session.flush()
            protein = session.query(Protein).filter(and_(Protein.sequence_key==sequence.id,Protein.experiment_key==E_ID)).first()
            if protein==None:
                protein=Protein()
                protein.experiment_key=E_ID
                protein.protein_type="phylogeny"
                protein.sequence_key = sequence.id
                protein.probability = 0
                protein.comment = "auto added amnh families"
                protein.file_key = 0
                protein.parse_key = 0
                protein.gene_key = 0
                session.add(protein)
                session.flush()
                debug("Added",protein)
            if not family in sequence.families:
                sequence.families.append(family)
                session.flush()
    session.commit()        
    session.close()
    debug("closed")
    return None
    
def main(*args):
    # Clear the db mappings, to use one engine per subprocess
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS))
    runtime().debug("Using processor",pool)
    pool.make_tasks(lambda: args)
    consume(pool.run(_import_family))
    
def _do(argv):
    r = runtime()
    r.description("""
    template.py [-options] args
    template script.
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing"))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
