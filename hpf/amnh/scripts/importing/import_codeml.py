#!/usr/bin/env python

import sys, os
import getopt
from hpf.runtime import runtime, Option
from hpf.utilities import consume
from hpf.hddb.db import *
from itertools import imap
from hpf.amnh.codeml import PositiveSelectionParser
            
FORMAT = "format"
session = None
family = None
tree = None

def merge(selection):
    selection.family_key = family
    selection.tree_key = tree
#    try:
#        s.sequence_key = int(record.id)
#    except:
#        seq = session.query(Sequence).filter(Sequence.sha1==sha1(str(record.seq).replace("-","")).hexdigest()).first()
#        if seq:
#            s.sequence_key = seq.id
#        else:
#            runtime().debug("Can't find sequence_key for",record.id, str(record.seq).replace("-",""))
#            return
    selection = session.merge(selection)
    runtime().debug("Merged",selection)
    return selection 

def main(fam,t, *files):
    global family, tree, session
    session = Session()
    family = int(fam)
    tree = int(t)
    
    for file in files:
        #runtime().set_debug(1)
        runtime().debug("Using file",file)
        with open(file) as handle:
            ps = PositiveSelectionParser().parse(handle)
        count = consume(imap(merge,ps))
        runtime().debug("Found",count,"sites")
            
    session.commit()
    session.close()
    
    
def _do(argv):
    r = runtime()
    r.description("""
    import_codeml.py [-options] familykey treekey *codeml_output
    Import codeml positive selection output.
    """)
    #r.add_option(Option(FORMAT, "f", description="Alignment Format", default="fasta"))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
