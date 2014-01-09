#!/usr/bin/env python

import sys, os
import getopt
from hpf.runtime import runtime, Flag
from hpf.utilities import consume
from hpf.hddb.db import *
from Bio import SeqIO
from hashlib import sha1
from sqlalchemy.sql import and_, or_, select, func

session = None
experiment = None
fasta = None

def digest(record):
    return sha1(str(record.seq)).hexdigest()

def add(record):
    record._digest = digest(record)
    z = session.query(Sequence).filter(Sequence.sha1==record._digest).first()
    if not z:
        z = Sequence(sequence = str(record.seq),
                     sha1= record._digest)
        session.add(z)
        session.flush()
        runtime().debug("Added sequence", record.id, z.id)
    else:
        runtime().debug("Found sequence", record.id, z.id)
    record._hddb = z
    record = ac(record)
    record = protein(record)
    record.id = str(record._ac.protein_key)
    return record
    
def ac(record):
    acc = session.query(SequenceAc).filter(and_(SequenceAc.sequence_key==record._hddb.id,
                                               SequenceAc.db=='amnh',
                                               SequenceAc.ac==record.id)).first()
    if not acc:
        acc = SequenceAc(sequence_key=record._hddb.id,
                        gi = None,
                        db = 'amnh',
                        ac = record.id,
                        ac2 = None,
                        description = fasta,
                        taxonomy_id = experiment.taxonomy_id)
                        
        session.add(acc)
        session.flush()
        runtime().debug("Added ac", acc)
    else:
        runtime().debug("Found ac", acc)
    
    record._ac = acc
    return record

def protein(record):
    if record._ac.protein_key == None:
        protein = Protein(experiment_key=experiment.id,
                          protein_type='amnh',
                          sequence_key=record._hddb.id,
                          probability=0,
                          comment="AMNH protein families",
                          file_key=0,
                          parse_key=0,
                          gene_key=0,
                          insert_data=func.CURRENT_DATE()
                          )
        session.add(protein)
        session.flush()
        record._ac.protein_key = protein.id
#        record._ac = session.merge(record._ac)
#        record._ac.protein_key = protein.id
        runtime().debug("Modified",session.is_modified(record._ac))
        session.add(record._ac)
        record._protein = protein
        runtime().debug("Added protein",record._protein)
    else:
        record._protein = session.query(Protein).get(record._ac.protein_key)
        runtime().debug("Found protein key",record._protein)
    return record
    
def map(function, sequence):
    result = []
    all = len(sequence)
    i = 0
    for item in sequence:
        i+=1
        result.append(function(item))
        if i%100 == 0:
            print "Finished %i, remaining %i" % (i, all-i)
    return result
        
def main(file, experiment_id):
    global fasta, session, experiment
    
    #runtime().set_debug(1)
    runtime().debug("Using file",file)
    session = Session()
    fasta = os.path.basename(file)
    
    experiment = session.query(Experiment).get(experiment_id)
    with open(file) as handle:
        records = list(SeqIO.parse(handle,"fasta"))
    runtime().debug("Found %i records" % len(records))
    runtime().indent()
    records = map(add,records)
    session.close()
    runtime().unindent()
    runtime().debug("Writing %i records" % len(records))
    with open(file+".hpf","w") as handle:
        SeqIO.write(records, handle, "fasta")
    
    
def _do(argv):
    r = runtime()
    r.description("""
    import_proteome.py [-options] fasta experiment_id
    template script.
    """)
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
