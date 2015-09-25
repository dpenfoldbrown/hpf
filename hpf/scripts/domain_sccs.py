#!/usr/bin/env python
# Uses domain tables to generate SCOP domain predictions for pdbblast and fold_recognition hits.

import sys
import getopt
import MySQLdb
import string
import multiprocessing
from hpf.processing import MultiProcessor, MapProcessor
from hpf.pdb.scop import SCOPDomain
from hpf.hddb.db import DomainFoldableMap, Domain, Session
session = Session()


HELP = "help"
DEBUG = "debug"
SILENT = "silent"
UPLOAD = "upload"
PICKLE = "pickle"
EXPERIMENT = "experiment"
MCM = "mcm"
PDB = "pdb"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, PICKLE: None, EXPERIMENT:None, PDB:False, MCM:False}
DOMAIN_SCCS_TABLE_NAME = "hpf.domain_sccs"
FOLD_COVERAGE_MIN = 0.75

db = None
# A synchronized dictionary for multiprocessing
manager = None
pdb = None
chain_lengths = None

def pdb_sccs(pdbid, chain_matched):
    """Request scop domains for a given pdb"""
    # Make sure it gets cached
    _pdb_sccs_cache(pdbid)
    return pdb[pdbid]
    
def pdb_domain(task):
    """Request a domain entry be mapped to a pdb superfamily.
    @param task: Tuple (domain_key, parent_key, domain_type, pdb_id, chain, p_start, p_stop, seq_len)
    @note: seq_len The length of the domain sequence. 
    """
    
    domain_key, parent_key, domain_type, pdb_id, chain, p_start, p_stop, seq_len = task
    
    scop_domains = pdb_sccs(pdb_id, chain)
    #print scop_domains
    # Require 50% overlap with scop domain, or 50% of our domain to be part of the scop
    matches = [domain for domain in scop_domains if domain.overlap(pdb_id, chain, p_start, p_stop) > domain.length()*0.5]#seq_len*0.5]
    if len(matches) > 0:
        sccs = ", ".join([match.sccs for match in matches])
    else:
        sccs = None
    size = max(len(matches),1)
    confidence=1.0
    #if sccs is not None:
    #    print pdb_id,domain_key,sccs
    #print sccs
    return (domain_key, parent_key, domain_type, sccs, confidence, pdb_id, chain, size)


def psiblast_fold():
    global db
    if db==None:
        db = _connection()
    cursor = db.cursor()
    try:
        print "PSI-Blast and Fold Recognition Query"
        query = """
            select distinct d.domain_sequence_key, d.parent_sequence_key,
                d.domain_type, i.pdbId, r.chain, n.parent_start, n.parent_stop,
                length(s.sequence)
            from hpf.experiment e
            join hpf.protein p
            on e.id=p.experiment_key
            join hpf.domain d
            on p.sequence_key=d.parent_sequence_key
            join hpf.sequence s
            on d.domain_sequence_key=s.id
            join hpf.pdbSeqRes r
            on substring(d.parent_id FROM 4)=r.sequence_key
            join hpf.domainRegion n
            on d.id=n.domain_key
            join hpf.pdbIndex i
            on r.pdb_key=i.id
            where d.domain_type in ('fold_recognition','psiblast') 
        """

        if opts[EXPERIMENT]:
            query = query+" and e.id in (%s)" % opts[EXPERIMENT]
        print query
        cursor.execute(query)
        tasks = [(domain_key, parent_key, domain_type, pdb_id, chain, int(p_start), int(p_stop), seq_len) for domain_key, parent_key, domain_type, pdb_id, chain, p_start, p_stop, seq_len in cursor.fetchall()]
        keys = {}
        # Cartesian product returns many chains/pdbs for the same sequence
        # Filter these to one pdb/chain per domain by hashing
        for task in tasks:
            domain_key, parent_key, domain_type, pdb_id, chain, p_start, p_stop, seq_len = task
            keys[(domain_key,parent_key)] = task    
        tasks = keys.values()
    finally:
        cursor.close(); db.close(); db=None; cursor=None
    
    global pdb, manager, chain_lengths
    #manager = multiprocessing.Manager()
    #pdb = manager.dict()
    #chain_len = manager.dict()
    pool = MultiProcessor(raise_errors=False, modulus=100,processors=6)
    print "PSI-Blast and Fold Recognition Process"
    results = []
    for result in pool.run(pdb_domain, tasks):
        if isinstance(result, Exception):
            print result
        else:
            results.append(result)
    return results

def testCoverage(mcm_result):
	print mcm_result
	coverage_entries = session.query(DomainFoldableMap.fold_coverage).filter(DomainFoldableMap.domain_sequence_key == mcm_result[0]).filter(DomainFoldableMap.outfile_key == Domain.outfile_key).filter(Domain.domain_sequence_key == mcm_result[0]).filter(Domain.parent_sequence_key == mcm_result[1]).all()
	for fold_coverage in coverage_entries:
		print "coverage: ",fold_coverage
		if float(fold_coverage[0]) <= FOLD_COVERAGE_MIN:
			print "remove: ",fold_coverage
			return False
		else:
			print "keep: ",fold_coverage
			return True

	return True
				

def mcm():
    global db
    if db==None:
        db = _connection()
    cursor = db.cursor()
    try:
        print "MCM Query"
	#kdrew: this originally joined on d.domain_sequence_key=m.sequence_key but now is by outfile key
        query = """select distinct d.domain_sequence_key, d.parent_sequence_key, d.domain_type, m.sccs, m.probability from hpf.experiment e join hpf.protein p join hpf.domain d join hpf.mcm m 
            on e.id=p.experiment_key and p.sequence_key=d.parent_sequence_key and d.outfile_key = m.outfile_key 
            where d.domain_type in ('msa','pfam','unassigned','user_defined')"""
        if opts[EXPERIMENT]:
            query = query+" and e.id in (%s)" % opts[EXPERIMENT]
        cursor.execute(query)
        results = [(domain_key, parent_key, domain_type, sccs, confidence, None, None, 1) for domain_key, parent_key, domain_type, sccs, confidence in cursor.fetchall()]
	print results
	#kdrew: filters out results for domains that do not match the foldable sequence within a certain coverage percent
	results[:] = [result for result in results if testCoverage(result)]
	
    finally:
        cursor.close(); db.close(); db=None; cursor=None
    return results

def sort(results):
    """Picks highest confidence score for each domain"""
    #print "Sorting",len(results)
    all = []
    r = {}
    for result in results:
        if not result:
            raise Exception("None Result",result)
        domain_key, parent_key, domain_type, sccs, confidence, pdb_id, chain, size = result
        if r.has_key(domain_key):
            other = r[domain_key]
            if not size>other[7] and not confidence>other[4]:
                # Skip bad ones
                continue
        r[domain_key] = result
        
    all = [r[domain] for domain in r]
    #print "Domains",len(all)
    return all

def upload(results):
    
    def convert(results):
        for result in results:
            domain_key, parent_key, domain_type, sccs, confidence, pdb_id, chain, size = result
            #print domain_key, domain_type, size, sccs, confidence
            yield (domain_key, parent_key, domain_type, sccs, confidence, pdb_id, chain, size)
    global db
    if db==None:
        db = _connection()
    cursor = db.cursor()
    print "Uploading", len(results)
    print "SCCS for ", len([r for r in results if r[3] != None])
    cursor.executemany(
        """INSERT INTO """+DOMAIN_SCCS_TABLE_NAME+""" (domain_sequence_key, parent_sequence_key, domain_type, sccs, confidence, pdbid, chain, size)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s) ON DUPLICATE KEY UPDATE parent_sequence_key=VALUES(parent_sequence_key),sccs=VALUES(sccs),domain_type=VALUES(domain_type),confidence=VALUES(confidence),pdbid=VALUES(pdbid),chain=VALUES(chain),size=VALUES(size)""",
        convert(results))
    cursor.close(); db.close(); cursor=None; db=None;
        

def main():
    import os
    import pickle
    if opts[MCM]:
        m = mcm()
    else:
        m = []
    if opts[PDB]:
        p = psiblast_fold() 
    else:
        p = []    
    
    results = sort( m + p)#+sort(mcm())
    if opts[UPLOAD]:
        upload(results)   
    print "Finished"
    
def _pdb_sccs_cache(pdbid):
    """Looks up and caches a dictionary of chains and xranges for a pdbid"""
    global db,pdb,chain_lengths
    
    # Each process has to keep it's own global dictionary because of
    # concurrency issues.
    if pdb == None:
        pdb = {}
    if chain_lengths == None:
        chain_lengths = {}
    
    # Search for the pdb
    if not pdb.has_key(pdbid):
        if db==None:
            db = _connection()
        cursor = db.cursor()
        try:
            query = "select part_text, sccs from scop_cla where pdb='%s'" % pdbid
            cursor.execute(query)
            all = cursor.fetchall()
            # Parse parts and chains
            domains = []
            for txt, sccs in all:
                sccs = sccs.strip()
                if len(sccs) == 0:
                    sccs = None
                s = SCOPDomain(pdbid, sccs, part_text=txt)
                for part in s.parts:
                    # Determine the length of the chain for parts that might not
                    # have a start/stop defined.
                    if not chain_lengths.has_key((pdbid,part.chain)):
                        query = "select length(sequence) from pdbIndex i join pdbSeqRes r on i.id=r.pdb_key join sequence s on r.sequence_key=s.id where i.pdbId='%s' and r.chain='%s'" % (pdbid,part.chain)
                        #print query
                        cursor.execute(query)
                        chain_len = cursor.fetchone()[0]
                        chain_lengths[(pdbid,part.chain)]=chain_len
                    else:
                        chain_len = chain_lengths[(pdbid,part.chain)]
                    # Tell the par the length of the chains
                    part.chain_len = chain_len
                domains.append(s)
            pdb[pdbid]=list(set(domains))
        finally:
            cursor.close()
        


def _connection():
    #return MySQLdb.connect(host="127.0.0.1", user="pwinters", passwd="bonneaulab", db="hpf", port=13307)    
    #return MySQLdb.connect(host="localhost", user="kdrew", passwd="kdrew_nyu", db="hpf")
    return MySQLdb.connect(host="localhost", user="dpb", passwd="dpb_nyu", db="hpf")    

def _do(argv):
    try:
        opts, args = _options(argv)
    except Exception, e:
        print e
        _usage()
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "?dxqe:pm")
    global opts
    for o,a in _opts:
        if o in ('-?',HELP):
            _usage()
            sys.exit()
        if o == '-d':
            opts[DEBUG] = True
        if o == '-x':
            opts[UPLOAD] = False
        if o == '-q':
            opts[SILENT] = True
        if o == '-e':
            opts[EXPERIMENT] = a
        if o == '-m':
            opts[MCM] = True
        if o == '-p':
            opts[PDB]= True
    if not opts[MCM] and not opts[PDB]:
        opts[MCM] = True
        opts[PDB] = True
    return opts, args
    
def _usage():
    print "$ domain_sccs.py [-options]"
    print "Calculate SCOP classifications from GINZU psiblast and" 
    print "fold_recognition domains and from MCM scores."
    print "Updates domain_sccs table when re-run."
    print ""
    print "Options"
    print "-e experiment number"
    print "-p run psiblast/foldrecognition"
    print "-m run mcmdata"

    sys.exit(2)

if __name__=="__main__":
    _do(sys.argv[1:])
