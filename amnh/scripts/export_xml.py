#!/usr/bin/env python

import sys
import getopt
from hpf.runtime import runtime, Flag
from hpf.utilities import consume
from hpf.processing import processor, SYNCHRONOUS

from hpf.amnh.xml.structure import *
from hpf.xml import DefaultXMLGenerator

GZIP = "gzip"

def xml(id):
    from hpf.hddb.db import Session, Family
    session = Session()
    family = session.query(Family).get(id)
    filename = "%i.xml" % family.id

    if runtime().opt(GZIP):
        import gzip
        filename = "%s.gz" % filename
        handle = gzip.open(filename,"w")
    else:
        handle = open(filename,"w")

    try:
        doc = FamilyFeatureBuilder(
            lambda: DefaultXMLGenerator(handle,pretty=True),
            lambda handler: StructureFeatureProvider(handler),
            lambda handler: ColumnFeatureProvider(handler),
            lambda handler: IeaFeatureProvider(handler),
            lambda handler: SelectionFeatureProvider(handler)
            )
        doc.buildDocument(family)
    finally:
        handle.close()
    session.close()

def tasks():
    from hpf.hddb.db import Session, Family
    session = Session()
    ids = session.query(Family.id).filter(Family.manually_curated==0).all()
    session.close()
    return [i[0] for i in ids]

def main(*args):
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS))
    runtime().debug("Using processor",pool)
    pool.make_tasks(tasks)
    consume(pool.run(xml))
    
def _do(argv):
    r = runtime()
    r.description("""
    export_xml.py [-options] args
    Exports full structure/family reports for families to xml files. 
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing"))
    r.add_option(Flag(GZIP, "z", description="GZip the resulting files."))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
