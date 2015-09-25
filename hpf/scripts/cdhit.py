#!/usr/bin/env python

import sys
import getopt
from hpf.cdhit import CDHitParser, CREATE_TABLE, INSERT_ENTRY, insert_format
import MySQLdb

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
TABLE = "table"
EXPERIMENT = "experiment"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, TABLE: None, EXPERIMENT: None}

def main(file):

    def entries():
        parser = CDHitParser(file)
        for cluster in parser:
            _debug( "Cluster",cluster.id)
            for entry in cluster:
                _debug("\t",entry.id,entry.perc,entry.repr)
                yield entry

    if opts[UPLOAD]:
        with _connection() as cursor:
            upload(entries(),opts[TABLE],cursor,experiment=opts[EXPERIMENT],foreign_key=lambda entry: entry.id)
    else:
        for entry in entries():
            pass

def upload(entries,table,cursor,experiment=None,foreign_key=lambda entry: None):
    query = CREATE_TABLE % ("if not exists "+table)
    _debug(query)
    cursor.execute(query)
    query = "insert into %s %s" % (table,INSERT_ENTRY)
    def format(es):
        for entry in es:
            yield insert_format(entry,foreign_key=foreign_key(entry),experiment=experiment)
    for e in format(entries):
        cursor.executemany(query,[e])

def _connection():
    return MySQLdb.connect(user="patrick",passwd="patrick_nyu")

def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "dxqe:t:", [DEBUG,SILENT,UPLOAD,EXPERIMENT,TABLE])
    global opts
    for o,a in _opts:
        if o in ('-d',DEBUG):
            opts[DEBUG] = True
        if o in ('-x', UPLOAD):
            opts[UPLOAD] = False
        if o in ('-q',SILENT):
            opts[SILENT] = True
        if o in ('-e',EXPERIMENT):
            opts[EXPERIMENT] = a
        if o in ('-t',TABLE):
            opts[TABLE] = a

    if opts[UPLOAD]:
        assert opts[TABLE]!=None, "-t TABLE must be defined if uploading results"
    assert(len(args) == 1)
    return opts, args

def _debug(*args):
    if opts[DEBUG]:
        _print(*args)

def _print(*args):
    if not opts[SILENT]:
        print " ".join([str(a) for a in args])

def _usage():
    print "$ cdhit.py [-options] file"
    print "Load a cdhit file"
    print ""
    print "Options"


if __name__=="__main__":
    _do(sys.argv[1:])
