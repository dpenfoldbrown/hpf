#!/usr/bin/env python

import sys, os
import getopt
from hpf.runtime import runtime, Option
from hpf.utilities import consume
from hpf.hddb.db import *
from Bio import AlignIO
from hashlib import sha1

FORMAT = "format"
session = None
family = None
fasta = None

def merge(record):
    s = FamilySequence()
    s.family_key = family
    s.sequence_key = int(record.id)
#    try:
#        s.sequence_key = int(record.id)
#    except:
#        seq = session.query(Sequence).filter(Sequence.sha1==sha1(str(record.seq).replace("-","")).hexdigest()).first()
#        if seq:
#            s.sequence_key = seq.id
#        else:
#            runtime().debug("Can't find sequence_key for",record.id, str(record.seq).replace("-",""))
#            return
    s.alignment = str(record.seq)
    s.seed = True
    s = session.merge(s)
    runtime().debug("Merged",s.sequence_key,s.family_key, s.alignment)
    return s 
        
def main(file, fam):
    global family, session
    session = Session()
#    if isinstance(fam, basestring):
#        f = Family()
#        f.name = fam
#        f.experiment_key = 0
#        session.merge(f)
#        session.flush()
#        family = f.id
#    else:
    family = int(fam)
    
    #runtime().set_debug(1)
    runtime().debug("Using file",file)
    with open(file) as handle:
        from hpf.amnh.oid import phylip
        records = phylip(handle)._records
    
    runtime().debug("Found %i records" % len(records))
    from hpf.amnh.oid import index
    index(records)
    records = map(merge,records)
    session.commit()
    session.close()
    
    
def _do(argv):
    r = runtime()
    r.description("""
    import_alignment.py [-options] fasta family
    Import an alignment and attach to a family.
    """)
    r.add_option(Option(FORMAT, "f", description="Alignment Format", default="fasta"))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
