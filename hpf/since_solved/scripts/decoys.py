#!/usr/bin/env python
# Export all the decoys/astrals grouped by cluster id
import os
import sys
import getopt
import shutil
import MySQLdb
from hpf.pdb.astral import AstralList

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
DIRECTORY = "directory"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, DIRECTORY: os.getcwd()}

__hpf,__bddb = (None,None)

def export(cluster_id,sequence_key,structure_key,probability,astral_id,ini_rms):
    dest = os.path.join(opts[DIRECTORY],str(cluster_id))
    try:
        os.makedirs(dest)
    except:
        pass
    astral_file = AstralList().retrieve_astral_file(astral_id)
    shutil.copy(astral_file,dest)
    decoy_file = os.path.join(dest,"%i_%i.decoy" % (sequence_key,structure_key))
    rms_file = os.path.join(dest,"%i_%i.rms.%.1f" % (sequence_key,structure_key,ini_rms))
    mcm_file = os.path.join(dest,"%i_%i.mcm.%.3f" % (sequence_key,structure_key,probability))
    for f in [rms_file,mcm_file]:
        with open(f,"w") as stat:
            pass
    db = _bddb()
    cursor = db.cursor()
    cursor.execute("select id,uncompress(compress_file_content) from ddbCommon.structure where id=%i" % structure_key)
    key,content = cursor.fetchone()
    cursor.close()
    with open(decoy_file,"w") as decoy:
        print >>decoy, content
        decoy.flush()

def main():  
    query = "select a.cluster_id,m.sequence_key,m.structure_key,m.probability,b.subject,o.ini_rms from (select c.cluster_id,c.identifier as sequence_key,min(b.expect) as expect from since_solved s join cdhit_clstr c on s.cluster_id=c.cluster_id join blast b on c.identifier=b.query group by c.identifier) a join hpf.mcmData m on a.sequence_key=m.sequence_key join blast b on b.query=m.sequence_key and a.expect=b.expect join mammoth o on m.structure_key=o.prediction group by m.structure_key"
    db = _hpf(db="since_solved")
    cursor = db.cursor()
    _debug(query)
    cursor.execute(query)
    for cluster_id,sequence_key,structure_key,probability,astral_id,ini_rms in cursor.fetchall():
        export(cluster_id,sequence_key,structure_key,probability,astral_id,ini_rms)

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

def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "?dxqo:", ["help",DEBUG,SILENT,UPLOAD,DIRECTORY])
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
        if o in ('-o',DIRECTORY):
            opts[DIRECTORY] = a
    return opts, args

def _debug(*args):
    if opts[DEBUG]:
        _print(*args)

def _print(*args):
    if not opts[SILENT]:
        print " ".join([str(a) for a in args])

def _usage():
    print "$ decoys.py [-options] args"
    print "Output decoys by cluster, sequence, structure, and rmsd with their astral counterparts."
    print ""
    print "Options"
    print "-o output directory"
    print ""


if __name__=="__main__":
    _do(sys.argv[1:])
