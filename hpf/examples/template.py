#!/usr/bin/env python

import sys
import getopt
from hpf.runtime import runtime, Flag
from hpf.utilities import consume
from hpf.processing import processor, SYNCHRONOUS

def main(*args):
    # calls processor with a the keyword argument 'synchronous'. runtime() creates a new Runtime obj. or accesses
    # the existing runtime object, and .opt(SYNCHRONOUS) returns the value of the synchronous option of the Runtime
    # object. 
    # processor returns a Map, SGEArray, or PBSArrayProcessor object. NOTE: currently, PBSArrayProcessor is not fully implemented.
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS))
    runtime().debug("Using processor",pool)
    pool.make_tasks(None)
    consume(pool.run(None))
    
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
