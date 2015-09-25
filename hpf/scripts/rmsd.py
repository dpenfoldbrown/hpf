#!/usr/bin/env python

import sys, os
import getopt
import hpf
import hpf.pdb
import hpf.processing
import MySQLdb
import Bio.PDB
import tempfile
from Bio.PDB.PDBIO import PDBIO
from hpf.pdb import Mammoth

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
SEQUENCES = "sequences"
EXPERIMENT = "experiment"
DECOYS = "decoys" 
opts = {DEBUG: False, SILENT: False, UPLOAD: True, SEQUENCES:None, EXPERIMENT:None, DECOYS:os.getcwd()}

_db=None

def rmsd(pdbid,chain,decoy_file,astral=False):
    """Calculate the RMS between the two structures based on sequence alignment."""
    if astral:
        pdb = hpf.pdb.get_astral(pdbid)
    else:
        pdb = hpf.pdb.get_pdb(pdbid)
    
    temp = tempfile.NamedTemporaryFile("w")
    try:
        # Save the pdbchain
        io = PDBIO()
        io.set_structure(pdb)
        io.save(temp, select=hpf.pdb.PDBChainSelector(chain))
        # Perform Mammoth alignment
        decoy_id = os.path.basename(decoy_file)
        cl = Mammoth.MammothCL(temp.name, decoy_file)
        mm = Mammoth.do_alignment(cl)
        return mm 
    finally:
        temp.close()

def do_it(arg):
    pred_code,seq_key,cluster_index,probability,pdbid,chain = arg
    decoy_path = os.path.join(opts[DECOYS],pred_code[0:2],pred_code,"decoy_%s.pdb"%cluster_index)
    mm = rmsd(pdbid,chain,decoy_path)
    return (seq_key,pred_code,cluster_index,probability,pdbid,chain,mm)


def tasks():
    con = _connection()
    query = """
select f.prediction_code,m.sequence_key,m.cluster_center_index,m.probability,q.pdbid, q.chain
from hpf2.gos g join
(select max(m.probability) as probability, g1.cluster_id, d.pdbid, d.chain from 
hpf.mcmData m join hpf2.gos g1 on m.sequence_key=g1.domain_sequence_key 
join hpf2.gos g2 on g1.cluster_id=g2.cluster_id
join hpf.domain_sccs d on g2.domain_sequence_key=d.domain_sequence_key
where d.domain_type in ("fold_recognition","psiblast")
group by g1.cluster_id order by probability desc ) q
on g.cluster_id=q.cluster_id
join hpf.mcmData m on m.sequence_key=g.domain_sequence_key and m.probability=q.probability
join hpf.filesystemOutfile f on m.sequence_key=f.sequence_key;
"""
    _debug(query)
    cursor = con.cursor()
    cursor.execute(query)
    for pred_code,seq_key,cluster_index,probability,pdbid,chain in cursor.fetchall():
            yield (pred_code,seq_key,cluster_index,probability,pdbid,chain)
    cursor.close(); con.close()

def calculate():
    t = list(tasks())
    pool = hpf.processing.MultiProcessor(7, modulus=1, raise_errors=True)
    for result in pool.run(do_it, t):
        seq_key,pred_code,cluster_index,probability,pdbid,chain,mm = result
        if not mm.none():
            yield result
    
def insert(results):
    con = _connection()
    cursor = con.cursor()
    upload = [(str(cluster_index), str(seq_key), str(pdbid), str(chain), str(mm.ini_psi), str(mm.ini_rms), str(mm.end_psi), str(mm.end_rms))for seq_key,pred_code,cluster_index,probability,pdbid,chain,mm in results]
    query = """INSERT IGNORE INTO hpf2.mammoth 
        (cluster_center_index, sequence_key, pdbid, chain, ini_psi, ini_rms, end_psi, end_rms)
        VALUES (%s,%s,'%s','%s',%s,%s,%s,%s)"""
    _debug(query)
    cursor.executemany(query, upload)

def display(results):
    for result in results:
        seq_key,pred_code,cluster_index,probability,pdbid,chain,mm = result
        print pred_code, cluster_index, seq_key, probability,pdbid,chain 
        print "\tini_psi",mm.ini_psi,"ini_rms",mm.ini_rms, "end_psi",mm.end_psi,"end_rms",mm.end_rms

def main():
    results = calculate()
    if opts[UPLOAD]:
        insert(results)
    else:
        display(results)

def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _connection(db=None):
    return MySQLdb.connect(read_default_file="~/.my.cnf")

def _options(argv):
    _opts, args = getopt.getopt(argv, "dxq", [DEBUG,SILENT,UPLOAD])
    global opts
    for o,a in _opts:
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
