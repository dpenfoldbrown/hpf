#!/usr/bin/env python

import sys
import os
import getopt
import subprocess

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "upload"
TASK_ID = "SGE_TASK_ID"
opts = {DEBUG: False, SILENT: False, UPLOAD: True}

def main(input_file):
    i = 0
    f = open(input_file)
    #ex >ddb017261348 unkn|jgi|Chlre3|179525|jgi|Chlre3|179525|jgi|Chlre3|179525
    try:
        for line in f:
            line=line.strip()
            if line.startswith(">"):
                i+=1
                if i==opts[TASK_ID]:
                    sequence_key = int(line.split()[0][4:])
                    break
    finally:
        f.close()
    print "Sequence_key: %i" % sequence_key
    shell = "./DDB/ddb.pl -mode sequence -submode process -sequence_key %i" % sequence_key
    print shell
    retcode = subprocess.call(shell, cwd=os.getcwd(), shell=True)
    sys.exit(retcode)
    

def _do(argv):
    try:
        opts, args = _options(argv)
    except Exception, e:
        print e
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "dxqt:")
    global opts
    if os.environ.has_key(TASK_ID):
        opts[TASK_ID]=int(os.environ[TASK_ID])
    for o,a in _opts:
        if o == '-d':
            opts[DEBUG] = True
        if o == '-x':
            opts[UPLOAD] = False
        if o == '-q':
            opts[SILENT] = True
        if o == '-t':
            opts[TASK_ID] = int(a)
    return opts, args

def _usage():
    print "$ ginzu.py [-options] args"
    print "Runs ginzu on the nth sequence."
    print ""
    print "Options"
    print "-t Task ID, Otherwise taken from SGE_TASK_ID"
    
if __name__=="__main__":
    _do(sys.argv[1:])
