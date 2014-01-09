#!/usr/bin/env python
# Perform rmsd analysis on since solved

import os
import sys
import getopt
import MySQLdb
from tempfile import NamedTemporaryFile
from hpf import processing
from hpf.pdb.mammoth import Mammoth,MammothCL
from hpf.pdb.astral import AstralList

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
opts = {DEBUG: False, SILENT: False, UPLOAD: True}

__hpf,__bdbb = (None,None)

def rmsd(experiment, prediction):
    cl = MammothCL(experiment, prediction)
    mm = Mammoth(cl).run()
    return mm

def _rmsd(args):
    structure_key,astral_id = args
    """Retrieves the decoy from Lars' database and the astral structure from the web."""
    _debug(astral_id,structure_key)
    astral = AstralList().retrieve_astral_file(astral_id)
    db = _bddb()
    cursor = db.cursor()
    cursor.execute("select id,uncompress(compress_file_content) from hpf.structure where id=%i" % structure_key)
    key,content = cursor.fetchone()
    cursor.close()
    with NamedTemporaryFile() as temp:
        print >>temp, content
        temp.flush()
        mm = rmsd(astral,temp.name)
        mm.experiment=astral_id
        mm.prediction=structure_key
    return mm

def _astral(astral_id):
    AstralList().retrieve_astral_file(astral_id)
    
def _tasks():
    db = _hpf(db="since_solved")
    cursor = db.cursor()
    query = "select m.structure_key,b.subject from (select c.identifier as sequence_key,min(b.expect) as expect from since_solved s join cdhit_clstr c on s.cluster_id=c.cluster_id join blast b on c.identifier=b.query group by c.identifier) a join since_solved.mcmdata_redux m on a.sequence_key=m.sequence_key join blast b on b.query=m.sequence_key and a.expect=b.expect group by m.structure_key"
    _debug(query)
    cursor.execute(query)
    tasks = cursor.fetchall()
    cursor.close()
    db.close()
    
    global __hpf,__bddb
    __hpf = None
    __bddb = None
    return tasks

def _upload(mm):
    db = _hpf(db="since_solved")
    cursor = db.cursor()
    query = "insert into mammoth_redux (prediction, experiment, ini_psi, ini_rms, end_psi, end_rms, zscore, evalue) values (%s,%s,%s,%s,%s,%s,%s,%s) on duplicate key update ini_psi=values(ini_psi),ini_rms=values(ini_rms),end_psi=values(end_psi),end_rms=values(end_rms),timestamp=NOW(),zscore=values(zscore),evalue=values(evalue)"
    cursor.executemany(query,[(mm.prediction,mm.experiment,mm.ini_psi,mm.ini_rms,mm.end_psi,mm.end_rms,mm.zscore,mm.evalue)])
    cursor.close()
    _debug(mm)

def _bddb():
    global __bddb
    if not __bddb:
        __bddb = MySQLdb.connect(host="127.0.0.1",user="pwinters",passwd="bonneaulab",port=13307,db="pwinters")
    return __bddb

def _hpf(db=None):
    global __hpf
    if not __hpf:
        __hpf = MySQLdb.connect(db=db,read_default_file="~/.my.cnf")
    return __hpf

def main():
    tasks = _tasks()
    _debug(len(tasks)," tasks")
    if opts[DEBUG]:
        for task in tasks:
            _upload(_rmsd(task))
    else:
        pool = processing.MultiProcessor(raise_errors=True,processors=8)
        for r in pool.run(_astral,list(set([astral_id for struct,astral_id in tasks]))):
            pass
        pool.reset()
        for r in pool.run(_rmsd,tasks,_upload):
            # Consume the generator
            pass

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
