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

from hpf.processing import array_processor
from hpf import hddb

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
EXPERIMENT="experiment"
TASK_PICKLE = "task_pickle"
PBS = "pbs"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, EXPERIMENT:None,PBS:False,
        TASK_PICKLE:"/home/patrick/export_ginzu.tasks"}

_db = None

def _tasks():
    global _db
    db = _connection()
    cursor = db.cursor()
    query = """
        SELECT DISTINCT p.sequence_key FROM hpf.experiment e 
        join hpf.protein p on e.id=p.experiment_key
        join hpf.domain d on p.sequence_key=d.parent_sequence_key
        join hpf.yeastrcu y on p.sequence_key=y.sequence_key 
        left outer join hpf._yrc_ginzu s on y.yrc_sequence_key=s.yrc_sequence_key 
        where s.yrc_sequence_key is NULL"""
    #if opts[EXPERIMENT]:
    #    query += " where p.experiment_key=%i" % opts[EXPERIMENT]
    _debug(query)
    cursor.execute(query)
    tasks = [result[0] for result in cursor.fetchall()]
    _print("Loaded ",len(tasks)," sequences")
    cursor.close(); cursor = None; _db.close(); _db = None
    return tasks

def protein(sequence_key):
    #subprocess.check_call("nohup `which Xvfb` :20 -screen 0 1024x768x24 &> /dev/null & sleep 7",shell=True)
    print >>sys.stderr, os.getenv("HOSTNAME")+" "+os.getenv("SGE_TASK_ID")
    hddb.tunnel()
    #-Djava.awt.headless=true
    #DISPLAY=localhost:20.0
    cmd = "ddb.pl -mode export -submode ginzu2yeastrc -ginzu 1 -parent_sequence_key %i" % sequence_key
    subprocess.check_call(cmd,shell=True)
    return sequence_key

def main():
    processor = array_processor(opts[TASK_PICKLE],pbs=opts[PBS])
    processor.make_tasks(_tasks)
    for r in processor.run(protein):
        pass

def _connection(db="hpf"):
    global _db
    if _db == None:
        _db = MySQLdb.connect(host="127.0.0.1",port=13307,user="pwinters",passwd="bonneaulab")
    return _db

def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "?dxqe:p", [HELP,DEBUG,SILENT,UPLOAD,EXPERIMENT,PBS])
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
        if o in ('-p',PBS):
            opts[PBS] = True
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
