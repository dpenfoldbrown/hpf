#!/usr/bin/env python

import sys, os
import getopt
import MySQLdb

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "upload"
EXPERIMENT = "experiment"
SEQUENCES = "sequences"
CL_SEQUENCE = "cl_sequence"
NUMBER = "number"
OUTPUT = "output"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, EXPERIMENT: None, SEQUENCES: None, CL_SEQUENCE: None, NUMBER: 1, OUTPUT:os.getcwd()}

_db = None

def sequences(seqs=None):
    _print("Loading Sequences")
    if seqs == None:
        seqs = []
    db = _connection()
    cursor = db.cursor()
    if opts[SEQUENCES]:
        f = open(opts[SEQUENCES])
        ids = [line.strip() for line in f if line.strip() != "" and not line.startswith("#")]
        cursor.execute("select s.id, s.sequence from hpf.sequence s where s.id in (%s)" % ",".join(ids))
    elif opts[CL_SEQUENCE]:
        cursor.execute("select s.id, s.sequence from hpf.sequence s where s.id = %s" % opts[CL_SEQUENCE])
    elif opts[EXPERIMENT]:
        _debug("Experiment: ",opts[EXPERIMENT])
        cursor.execute("select s.id, s.sequence from hpf.protein p join hpf.domain d join hpf.sequence s on p.sequence_key=d.parent_sequence_key and d.domain_sequence_key=s.id where p.experiment_key=%i" % opts[EXPERIMENT])
    
    seqs += [(int(id),sequence) for id,sequence in cursor.fetchall()]
    cursor.close()
    _debug("Loaded %i sequences" % len(seqs))
    return seqs

def mcm(seqs, number):
    _print("Loading Best %i MCM scores" % number)
    db = _connection()
    cursor = db.cursor()
    query = "select m.probability, m.sequence_key, m.cluster_center_index, m.structure_key from hpf.mcmData m where m.sequence_key in (%s) order by m.sequence_key" % ",".join([str(id) for id, seq in seqs])
    _debug(query)
    cursor.execute(query)
    
    mcms = {}    
    for prob, key, cluster, struct in cursor.fetchall():
        if not mcms.has_key(key):
            mcms[key] = []
        mcms[key].append((prob,cluster, struct))
    cursor.close()
    
    results = []
    for key in mcms:
        mcms[key].sort(reverse=True)
        results += [(key,c,s) for p,c,s in mcms[key][:number]]
        
    return results
    
def decoys(mcms, outdir):
    _print("Writing %i Decoys" % len(mcms))
    db = _connection()
    cursor = db.cursor()
    ret_dict = {}
    for seq,cluster,struct in mcms:
        cursor.execute("select uncompress(compress_file_content),1 from hpf.structure where id=%i" % struct)
        dir = os.path.join(outdir,str(seq))
        if not os.path.exists(dir):
            os.mkdir(dir)
        f = open(os.path.join(dir,"decoy_%i.pdb" % struct), "w")
        content,blah = cursor.fetchone()
        _debug(content)
        print >>f, content
	ret_dict[struct] = f.name
    return ret_dict

def main(*args):
    #seqs = sequences(args)
    seqs = sequences()
    mcms = mcm(seqs, opts[NUMBER])
    decoys(mcms, opts[OUTPUT])
        
    
def _connection():
    global _db
    if _db == None:
        _db = MySQLdb.connect(host="127.0.0.1", user="pwinters", passwd="bonneaulab", db="hpf", port=3307)
    return _db

def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "dxqe:s:S:o:n:")
    global opts
    for o,a in _opts:
        if o == '-d':
            opts[DEBUG] = True
        if o == '-x':
            opts[UPLOAD] = False
        if o == '-q':
            opts[SILENT] = True
        if o == '-e':
            opts[EXPERIMENT] = int(a)
        if o == '-s':
            assert(os.path.exists(a))
            opts[SEQUENCES] = a
        if o == '-S':
            opts[CL_SEQUENCE] = a
        if o == '-o':
            assert(os.path.exists(a))
            opts[OUTPUT] = a
        if o == '-n':
            opts[NUMBER] = int(a)
    return opts, args

def _debug(*args):
    if opts[DEBUG]:
        _print(*args)

def _print(*args):
    if not opts[SILENT]:
        print " ".join([str(a) for a in args])

def _usage():
    print "$ extract_decoys.py [-options] *sequence_keys"
    print "Extract decoys for a given set of sequences, a sequence"
    print "    file, or experiment number."
    print ""
    print "Options"
    print "-e experiment number to grab sequences for."
    print "-s file to load sequence keys from."
    print "-n number of best decoys to include."
    print ""
    print "Output"
    print "Exports decoys in the following format"
    print "    ./$SEQUENCE_KEY/decoy_$STRUCTURE_KEY.ent"


if __name__=="__main__":
    _do(sys.argv[1:])
