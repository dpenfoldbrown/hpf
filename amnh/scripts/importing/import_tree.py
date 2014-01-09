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

session = None

def add(file, family_key):
    with open(file) as handle:
        from Bio.Nexus.Trees import Tree as NexusTree
        tree = NexusTree(handle.read())
    dbtree = Tree()
    dbtree.family_key = family_key
    dbtree.text = tree.to_string(plain=False,plain_newick=True)
    dbtree.filename = os.path.basename(file)
    session.add(dbtree)
    session.flush()
    tree.name = dbtree.id
    for node in TreeNodeFactory().create(tree):
        session.add(node)

    
def main(file, family_key):
    # Clear the db mappings, to use one engine per subprocess
    global session
    session = Session()
    add(file, family_key)
    session.commit()
    session.close()
    
def _do(argv):
    r = runtime()
    r.description("""
    import_tree.py [-options] tree_file family_key
    Import a tree into the database.
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing"))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
