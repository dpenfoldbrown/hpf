#!/usr/bin/env python
# Run the YRC function update method on individual proteins

import sys
import getopt
import subprocess
import MySQLdb
from hpf.processing import SGEArrayProcessor

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
opts = {DEBUG: False, SILENT: False, UPLOAD: True}

_db=None

def main():
    pool = SGEArrayProcessor(task_pickle="/home/patrick/yrc_functions.tasks")
    pool.make_tasks(_tasks)
    for r in pool.run(process):
        pass

def process(parent_sequence_key):
    # Force the tunnel to be open
    cmd = "/home/patrick/local/projects/ginzu/washington.sh"
    subprocess.call(cmd, shell=True)
    cmd = "ddb.pl -mode export -submode ginzu2yeastrc -functions 1 -parent_sequence_key %i" % parent_sequence_key
    _print(cmd)
    subprocess.check_call(cmd,shell=True)

def _tasks():
    """Load the protein keys to process"""
    with _connection() as cursor:
        query = """
            select distinct sequence_key from bddb.go where evidence_code='KD'
            """
        _debug(query)
        cursor.execute(query)
        return [t[0] for t in cursor.fetchall()]

def _connection(db="bddb"):
    global _db
    if _db == None:
        _db = MySQLdb.connect(host="127.0.0.1",port=13307,db=db,user="pwinters",passwd="bonneaulab")
    return _db

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
    print "$ yrc_functions.py.py [-options] args"
    print "Farm the yrc function updates to the grid.  Processes each protein"
    print "sequence key individually."
    print ""
    print "Options"


if __name__=="__main__":
    _do(sys.argv[1:])
