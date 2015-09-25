#!/usr/bin/env python

import os
import sys
import getopt
import subprocess
from cStringIO import StringIO
from Bio import SeqIO
from hpf.processing import SGEArrayProcessor
from hpf import hddb

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
FASTA = "fasta"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, FASTA:None}

def main(fasta_file):
    opts[FASTA] = fasta_file
    pool = SGEArrayProcessor(task_pickle="/home/patrick/import_nr.tasks")
    pool.make_tasks(lambda: None)
    for r in pool.run(process):
        pass

def process(n):
    fasta_file = opts[FASTA]
    fasta_file+=".%i" % n-1
    file = os.path.abspath(fasta_file)
    cmd = "ddb.pl -mode update -submode nr -filename %s" % file
    hddb.tunnel()
    _debug(cmd)
    subprocess.check_call(cmd,shell=True)
    cmd = "ddb.pl -mode update -submode nr"

def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "?dxq", [DEBUG,SILENT,UPLOAD])
    global opts
    for o,a in _opts:
        if o in ('-?',HELP):
            _usage()
            sys.exit()
        if o in ('-d',DEBUG):
            opts[DEBUG] = True
        if o in ('-x', UPLOAD):
            opts[UPLOAD] = False
        if o in ('-q',SILENT):
            opts[SILENT] = True
    return opts, args

def _debug(*args):
    if opts[DEBUG]:
        _print(*args)

def _print(*args):
    if not opts[SILENT]:
        print " ".join([str(a) for a in args])

def _usage():
    print "$ template.py [-options] args"
    print "Template for quick scripts."
    print ""
    print "Options"


if __name__=="__main__":
    _do(sys.argv[1:])
