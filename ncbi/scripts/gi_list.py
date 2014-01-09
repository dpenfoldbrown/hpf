#!/usr/bin/env python
import os
import sys
import getopt
import MySQLdb
from Bio import Entrez
from hpf.processing import SGEArrayProcessor
from hpf import ncbi


DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
EXPERIMENT = "experiment"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, EXPERIMENT:None}

_db = None

def main():
    # Override for single experiment
    if opts[EXPERIMENT]:
        id = _tasks()[0]
        return organism(id)
    
    pool = SGEArrayProcessor("/home/patrick/gi_list.tasks")
    pool.make_tasks(_tasks)
    for r in pool.run(organism):
        pass

def organism(taxonomy_id):
    file = "%i.gi" % taxonomy_id
    _debug(file)
    with open(file,"w") as out:
        for gi in ncbi.gis(taxonomy_id):
            _debug(gi)
            print >>out, "%i" % gi

def _tasks():
    with _connection() as cursor:
        query = "select distinct taxonomy_id from hpf.experiment where taxonomy_id is not NULL"
        if opts[EXPERIMENT]:
            query += " and id=%i" % opts[EXPERIMENT]
        _debug(query)
        cursor.execute(query)
        return [result[0] for result in cursor.fetchall()]

def _connection(db="hpf"):
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
    _opts, args = getopt.getopt(argv, "?dxqe:", [DEBUG,SILENT,UPLOAD])
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
        if o in ('-e',EXPERIMENT):
            opts[EXPERIMENT] = int(a)
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
