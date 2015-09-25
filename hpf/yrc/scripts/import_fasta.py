#!/usr/bin/env python

import os
import sys
import getopt
import subprocess
import MySQLdb
from cStringIO import StringIO
from Bio import SeqIO
from hpf.processing import SGEArrayProcessor
from hpf.hddb import tunnel

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
FASTA = "fasta"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, FASTA:None}
_db=None

def main(fasta_file):
    pool = SGEArrayProcessor()
    pool.make_tasks(lambda: None)
    for r in pool.run(process):
        pass

def _upload(seq_key,prot_keys,ddb_key):
    query = """
    insert ignore into hpf.yeastrcu 
    (yrc_sequence_key, yrc_protein_key, sequence_key) 
    values (%s,%s,%s)
    """
    _debug(query)
    db = _connection()
    cursor = db.cursor()
    try:
        for prot_key in prot_keys:
            values = (seq_key,prot_key,ddb_key)
            _debug(values)
            cursor.execute(query, values)
    finally:
        cursor.close()


def process(n):
    fasta_file = opts[FASTA]
    fasta_file+=".%i" % (int(n)-1)
    with open(fasta_file) as io:
        for record in SeqIO.parse(io,"fasta"):
            try:
                seq_key,prot_keys = keys(record)
                seq = record.seq.tostring()
                ddb_key = get_ddb_key(seq)
                _upload(seq_key,prot_keys,ddb_key)
            except Exception, e:
                _debug(e)
                continue
    
    
def get_ddb_key(seq):
    query = "select id from hpf.sequence where sha1=sha1('%s')" % seq.strip().replace(" ","").upper()
    with _connection() as cursor:
        _debug(query)
        cursor.execute(query)
        id = [r[0] for r in cursor.fetchall()]
        assert len(id)==1
        id = id[0]
        return id

def get_record(fasta_file, n):
    """
    Process the n-th record from the fasta file.
    @param n: Number record to process from fasta. Starts with 1.
    @deprecated: Not based on N records, but N fasta files
    """
    raise Exception("NO MORE")
    i = 0
    i = int(2*n)
    cmd = "head -%i %s" % (i,fasta_file)
    out,err = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).communicate()
    io = StringIO()
    #io.write(out)
    with open(fasta_file) as handle:
        for line in handle:
            i+=1
            if n*2==i or n*2-1==i:
                io.write(line)
            if i>n*2:
                break

    #cmd = "head -%i %s | tail -2" % (i,fasta_file)
    #out,err = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).communicate()
    #_debug(out)
    # Flush and point to beginning of handle for reading
    io.flush()
    io.reset()
    _debug(io.getvalue())
    records = list(SeqIO.parse(io, "fasta"))
    assert len(records)==1
    return records[0]
    

def keys(record):
    seq_key=record.name
    prot_keys = [int(s) for s in record.description.split()[1].split("|")]
    return (seq_key, prot_keys)

def _connection(db="hpf"):
    global _db
    if _db == None:
        tunnel(5)
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
    if len(argv)>0:
        opts[FASTA]=argv[-1]
    return opts, args

def _debug(*args):
    if opts[DEBUG]:
        _print(*args)

def _print(*args):
    if not opts[SILENT]:
        print " ".join([str(a) for a in args])

def _usage():
    print "$ import_fasta.py [-options] fasta_file"
    print "Parses a YRC fasta file and puts entries into yeastrcu database table."
    print ""
    print "Options"


if __name__=="__main__":
    _do(sys.argv[1:])
