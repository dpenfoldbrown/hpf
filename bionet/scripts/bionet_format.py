#!/usr/bin/env python
#SELECTS FROM HPF.GO

import os
import sys

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
EXPERIMENT = "experiment"
# Defines the LLR with which we expect
LLR_CONSTANT = "llr_constant"
MCM_CONSTANT = "mcm_constant"
TASK_PICKLE = "task_pickle"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, LLR_CONSTANT:6.23, MCM_CONSTANT:0.8, EXPERIMENT:None,
        TASK_PICKLE:"/home/patrick/bionet_format.tasks"}

_db = None

def protein(protein):
    d = domain(protein)
    f = function(protein)
    g = gi(protein)
    d.update(f)
    d.update(g)
    d["parent_sequence_key"] = protein[0]
    return d

def gi(protein_key):
    db = _connection()
    cursor = db.cursor()
    query = "select distinct gi from sequenceAc where sequence_key=%i" % protein_key
    _debug(query)
    cursor.execute(query)
    gis = cursor.fetchall()
    _debug(gis)
    return {"gi":gis}
    

def domain(protein_key):
    """
    Selects the domain, structure, and ginzu info for a protein_key
    """
    db = _connection()
    cursor = db.cursor(MySQLdb.cursors.DictCursor)
    query = "select distinct size,domain_type,confidence,sccs from domain_sccs where parent_sequence_key=%i" % protein_key
    #_debug(query)
    cursor.execute(query)
    rows = cursor.fetchall()
    cursor.close()
    size = sum([row["size"] for row in rows])
    pdb_size = sum([row["size"] for row in rows if row["domain_type"].strip() in ["psiblast","fold_recognition"] ])
    rosetta_size = sum([row["size"] for row in rows 
               if row["domain_type"].strip() in ["pfam","msa","unassigned","user_defined"] 
               and row['confidence'] >= opts[MCM_CONSTANT]])
    pdb = pdb_size > 0
    rosetta = rosetta_size > 0
    for row in rows:
            _debug(row["domain_type"].strip(), row["domain_type"].strip() in ["pfam","msa","unassigned","user_defined"], row["size"],row["confidence"])

    
    #_debug("pdb_size:",pdb_size,"rosetta_size:",rosetta_size)
    coverage = (float(pdb_size)+float(rosetta_size))/float(size)
    sccs = ", ".join([row["sccs"] for row in rows if row["confidence"] >= opts[MCM_CONSTANT] and row["sccs"] != None])
    
    ginzu_counts = defaultdict(lambda: 0)
    for row in rows:
        ginzu_counts[row["domain_type"]] += row["size"]
    ginzu = ", ".join(["%s:%i" % (domain_type, ginzu_counts[domain_type]) for domain_type in ["psiblast","fold_recognition","pfam","msa","unassigned","user_defined"] if ginzu_counts[domain_type]>0])
    quality=0
    if rosetta and not pdb:
        quality=1
    if pdb and rosetta:
        quality=2
    if pdb and not rosetta:
        quality=3
    return {"size":int(size),"coverage":coverage,"pdb":pdb,"rosetta":rosetta,"sccs":sccs,"ginzu":ginzu,"quality":quality}
    
def function(protein_key):
    """Selects the 'best' function prediction for a protein"""
    db = _connection()
    cursor = db.cursor()
    query = "select b.pls_llr,b.base_llr,b.name,b.mf_acc from functionPredictions.bayes_golite_062009_3 b where b.parent_sequence_key=%i" % protein_key
    _debug(query)
    cursor.execute(query)
    best_value = None
    best_name = None
    best_acc = None
    best_pls = None
    best_base = None
    # The method of selection is difficult to decide on.
    # Currently we're balancing specificity with probability.
    for pls_llr,base_llr,name,mf_acc in cursor.fetchall():
        scaled_pls_llr = float(pls_llr)-opts[LLR_CONSTANT]
        prob_c = math.exp(scaled_pls_llr)/(1+math.exp(scaled_pls_llr))
        value = -float(base_llr)*prob_c
        if best_value == None or value > best_value:
            best_value = value
            best_name = name
            best_acc = mf_acc
            best_pls = pls_llr
            best_base = base_llr
    cursor.close()
    return {"mf_value":best_value, "mf_name":best_name, "mf_acc":best_acc, "mf_pls_llr":best_pls, "mf_base_llr":best_base}

def upload(result):
    #print opts[UPLOAD]
    #print result
    query = """INSERT INTO hpf.bionet_builder (gi, parent_sequence_key, size, coverage, pdb, rosetta, quality, sccs, ginzu, mf_value, mf_name, mf_acc, mf_pls_llr, mf_base_llr)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) 
        ON DUPLICATE KEY UPDATE 
            gi=VALUES(gi),
            parent_sequence_key=VALUES(parent_sequence_key),
            size=VALUES(size),
            coverage=VALUES(coverage),
            pdb=VALUES(pdb),
            rosetta=VALUES(rosetta),
            quality=VALUES(quality),
            sccs=VALUES(sccs),
            ginzu=VALUES(ginzu),
            mf_value=VALUES(mf_value),
            mf_name=VALUES(mf_name),
            mf_acc=VALUES(mf_acc),
            mf_pls_llr=VALUES(mf_pls_llr),
            mf_base_llr=VALUES(mf_base_llr)"""
    #values = []
    values = [result[col] for col in ["parent_sequence_key", "size", "coverage", "pdb", "rosetta", "quality", "sccs", "ginzu", "mf_value", "mf_name", "mf_acc", "mf_pls_llr", "mf_base_llr"]]
#        v = result[col]
#        if isinstance(v, str):
#            v = "'%s'" % v
#        elif isinstance(v, float):
#            v = "%f" % v
#        elif isinstance(v, int):
#            v = "%i" % v
#        elif v == None:
#            v = ""
#        values.append(v)
#         
    if opts[UPLOAD]:
        _debug("GI",result["gi"])
        _debug(values)
        
        if result["gi"]==None or len(result["gi"]) == 0:
            result["gi"] = [tuple([None])]
            _debug(result["gi"])
        
        db = _connection()
        cursor = db.cursor()
        cursor.executemany(query, [(gi)+tuple(values) for gi in result["gi"]])
        _debug(query)
        cursor.close()
    return result

def _tasks():
    global _db
    db = _connection()
    cursor = db.cursor()
    query = "select distinct s.parent_sequence_key from protein p join domain_sccs s on p.sequence_key=s.parent_sequence_key left outer join bionet_builder b on s.parent_sequence_key=b.parent_sequence_key where b.id is NULL"
    #query = "select distinct s.parent_sequence_key from protein p join domain_sccs s on p.sequence_key=s.parent_sequence_key"
    if opts[EXPERIMENT]:
        query += " where p.experiment_key=%i" % opts[EXPERIMENT]
    _debug(query)
    cursor.execute(query)
    tasks = cursor.fetchall()
    _print("Loaded ",len(tasks)," sequences")
    cursor.close(); cursor = None; _db.close(); _db = None
    return tasks

def main():
    processor = SGEArrayProcessor(raise_errors=True, modulus=1000,task_pickle=opts[TASK_PICKLE])
    processor.make_tasks(_tasks)
    for r in processor.run(protein,result=upload):
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
