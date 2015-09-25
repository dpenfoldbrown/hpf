#!/usr/bin/env python
#SELECTS FROM HPF.GO

import os
import sys
import subprocess
import getopt
import math
import MySQLdb
try:
    from collections import defaultdict
except:
    from hpf.utilities import defaultdict

from hpf.processing import SGEArrayProcessor,MultiProcessor

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
TASK_PICKLE = "task_pickle"
opts = {DEBUG: False, SILENT: False, UPLOAD: True,
        TASK_PICKLE:"/home/patrick/export_functions.tasks"}

_db = None

def _tasks():
    global _db
    db = _connection()
    cursor = db.cursor()
    query = "select distinct s.parent_sequence_key from protein p join domain_sccs s on p.sequence_key=s.parent_sequence_key where sccs is not NULL"
    #if opts[EXPERIMENT]:
    #    query += " where p.experiment_key=%i" % opts[EXPERIMENT]
    _debug(query)
    cursor.execute(query)
    tasks = cursor.fetchall()
    _print("Loaded ",len(tasks)," sequences")
    cursor.close(); cursor = None; _db.close(); _db = None
    return tasks

def protein(sequence_key):
    cmd = "ddb.pl -mode export -submode ginzu2yeastrc -function 1 -parent_sequence_key %i" % sequence_key
    subprocess.check_call(cmd,shell=True)
    return sequence_key

def main():
    processor = SGEArrayProcessor(raise_errors=True, modulus=1000,task_pickle=opts[TASK_PICKLE])
    processor.make_tasks(_tasks)
    for r in processor.run(protein):
        pass

def _connection(db="hpf"):
    global _db
    if _db == None:
        _db = MySQLdb.connect(host="carl.bio.nyu.edu",db=db,read_default_file="~/.my.cnf")
    return _db

def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "?dxqe:", [DEBUG,SILENT,UPLOAD,EXPERIMENT])
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
        if o in ('-e', EXPERIMENT):
            opts[EXPERIMENT] = int(a)
    return opts, args

def _debug(*args):
    if opts[DEBUG]:
        _print(*args)

def _print(*args):
    if not opts[SILENT]:
        print " ".join([str(a) for a in args])

def _usage():
    print "$ bionet_format.py [-options] args"
    print "Processes bionetbuilder attributes. Selects pre-adjusted function pls from hpf.go."
    print ""
    print "Options"
    print "-e experiment number"

if __name__=="__main__":
    _do(sys.argv[1:])
