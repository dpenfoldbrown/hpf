#!/usr/bin/env python
import sys, os
import getopt
import MySQLdb
from MySQLdb.cursors import DictCursor
import cPickle
from datetime import datetime
from hpf.function.metric import *
from hpf.function.term import *
from hpf.function.domain import *
from hpf.function.bayes import *

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "upload"
EXPERIMENT = "experiment"
SEQUENCES = "sequences"
SEQUENCE_TABLE = "sequence_table"
DOMAIN_SCCS_TABLE = "domain_sccs_table"
MCM = "mcm"
FR = "fold_recognition"
PSI = "psiblast"
SHELVE = "SHELVE"
PROCESSORS = "processors"
OLDER = "olderthan"
FUNCTION_DATABASE = "fndb"
FUNCTION_SUFFIX = "suffix"
TERM_SUFFIX = "tsuffix"
UPLOAD_DATABASE = "updb"
UPLOAD_TABLE = "uptbl"
MYGO_DATABASE = "mygo"
COMBINE = "combine"
BEST = "best"
MCMDATA = "mcmdata"
STRUCTURE_QUERY = "structure"
opts = {DEBUG: False, 
        SILENT: False, 
        UPLOAD: True, 
        EXPERIMENT: None, 
        SEQUENCES: None,
        SEQUENCE_TABLE: None,
        DOMAIN_SCCS_TABLE: "hpf.domain_sccs",
        MCM: False,
        FR: False,
        PSI: False,
        SHELVE: None,
        PROCESSORS: 8,
        FUNCTION_DATABASE:"functionTables",
        FUNCTION_SUFFIX:"_goLite_062009_3",
        TERM_SUFFIX:"_golite_062009",
        UPLOAD_DATABASE:"functionPredictions",
        UPLOAD_TABLE:None,
        OLDER: "NOW()",
        MYGO_DATABASE:"mygoLite_062009",
        COMBINE: False,
        BEST: False,
        MCMDATA: False,
        STRUCTURE_QUERY: None}

_bfp = None
_db = None

def main():
    from hpf.processing import MultiProcessor
    global _bfp,_db
    
    if opts[SEQUENCES]:
        f = open(opts[SEQUENCES])
        opts[SEQUENCES] = [line.strip().split()[0] for line in f if line.strip() != "" and not line.startswith("#")]
    
    seqs = defaultdict(lambda: [])
    for parent_key,domain_key in sequences(sequence_keys=opts[SEQUENCES], experiment_keys=opts[EXPERIMENT]):
        seqs[parent_key].append(domain_key)
    # The tasks are grouped by parent_sequence_key
    tasks = [(pkey,seqs[pkey]) for pkey in seqs]
    _print(len(tasks)," proteins for processing")

    
    lr,mi,fterms = tuple(metrics())
    print len(fterms), " functions known"
    
    # Be careful to close and serialize the persistent dictionaries
    try:
        _bfp = BayesFunctionPredictionDB(lr,mi, fterms)
        # Open a multi-processor for performing predictions on multiple
        pool = MultiProcessor(processors=opts[PROCESSORS],raise_errors=True)
        # Consume the pool generator, ignoring the results
        for r in pool.run(_predict, tasks, result=_upload,batches=100000):
            pass
    finally:
        for metric in (lr,mi):
            if hasattr(metric, 'close'):
                metric.close()
        #lr.close()
        #mi.close()
    

# Load the MI and LogRatios and FunctionTerms
def metrics():
    global _db
    connection = _connection(db=opts[FUNCTION_DATABASE])
    try:
        for filename,tablename,constructor in [("LogRatios","log_ratio",LogRatios),
                                               ("MutualInformation","mutual_info",MutualInformation)]:
            metric = constructor()#filename=file)
            load=True
            # Will create shelves if option is set
            if opts[SHELVE]:
                file = os.path.join(opts[SHELVE],filename+opts[FUNCTION_SUFFIX]+".shelve")
                load = os.path.exists(file)
            if load:
                print "Loading metric from database"
                metric.load_metric(connection,tablename+opts[FUNCTION_SUFFIX],
                                   query="select m.acc, m.acc2, m.metric from %s as m join mygolite_062009.term t on m.acc=t.acc and t.term_type='molecular_function'")
            
            yield metric
    finally:
        connection.close()
        _db = None
    
    connection = _connection(db="mygolite_062009")
    try:
        query = "select distinct t.acc,t.name,t.term_type from %s.%s l join mygoLite_062009.term t on l.acc=t.acc and t.term_type='molecular_function' where l.acc2 is NULL;" % (opts[FUNCTION_DATABASE],"log_ratio"+opts[FUNCTION_SUFFIX])
        fterms = FuncTerms(connection,query=query)
        yield fterms
    finally:
        connection.close()
        _db = None 

def _predict(keys):
    """Process each protein and its domains separately."""
    global _bfp
    # grouped by parent_sequence_key, and then a list of domains
    parent_key, domain_keys = keys
    results = []
    # Create a domain from each parent_key,domain_key combo
#    for domain in domains([(parent_key, dkey) for dkey in domain_keys], fr=opts[FR], mcm=opts[MCM], psi=opts[PSI]):
    # temporary override for mcmdata
    if opts[MCMDATA]:
        structure_query="select m.experiment_sccs, m.probability, 'unassigned' as domain_type from mcmdata m join filesystemOutfile f on m.sequence_key=f.sequence_key where m.sequence_key=%i"
    elif opts[STRUCTURE_QUERY] != None:
        structure_query = opts[STRUCTURE_QUERY]
    else:
        structure_query="select sccs, confidence, domain_type from "+opts[DOMAIN_SCCS_TABLE]+" where domain_sequence_key=%i"
    for domain in domains([dkey for dkey in domain_keys], parent_key, fr=opts[FR], mcm=opts[MCM], psi=opts[PSI], structure_query=structure_query):#,:
        #print domain
        #print ['%s' % term.get_id() for term in domain.proc_terms] 
        #print ['%s' % term.get_id() for term in domain.loc_terms] 
        #print domain.structures
        results.append(_bfp.predictDomain(domain,lower_bound=-3))
    return results

def _upload(domains):
    if opts[UPLOAD] and domains != None and len(domains) > 0:
        connection = _connection()
        table = opts[UPLOAD_TABLE]
        for domain in domains:
            cursor = connection.cursor()
            #print "RESULTS",domain.predictions
            _bfp.uploadPredictions(cursor, domain, table)
            cursor.close()
    return None

def domains(domain_keys,parent_key, fr=False, mcm=False, psi=False, structure_query="select sccs, confidence, domain_type from "+opts[DOMAIN_SCCS_TABLE]+" where domain_sequence_key=%i"):
    """Loads domains and structures based on types requested.
    This will combine superfamilies or treat them as individual domains
    according to options."""
    term_table="hddb_IEA.hddb_IEA"+opts[TERM_SUFFIX]
    assert(fr | mcm | psi)
    f,m,p = (FRDomain, MCMDomain, PSIDomain)
    keep = {f:fr,m:mcm,p:psi}
    ginzu = {"fold_recognition":f, "msa":m, "pfam":m, "psiblast":p, "unassigned":m, "user_defined":m}

    
    # Loops over domain_keys... parent_key always the same.
    for domain_key in domain_keys:
        # Get all structures
        structs = sorted([Structure(sf,confidence,domain_type=domain_type) for sf,confidence,domain_type in _structures(domain_key,query=structure_query,fr=fr, mcm=mcm, psi=psi)],reverse=True)
        if opts[COMBINE]:
            # Keep superfamilies as a single set, but combine them
            set = Structures(structs)
            set.combine_multiple_sf()
            set = [set]
        elif opts[BEST]:
            # Pick the top score from the set
            set = [Structures([structs[0]])] if len(structs) >0 else []
        else:
            # Split them up into separate sets
            set = [Structures([s]) for s in structs]
        db = _connection()
        for structs in set:
            if len(structs)==0:
                continue
            else:
                domain_type=structs[0].domain_type
            domain = ginzu[domain_type](psk=parent_key, dsk=domain_key, evidence_codes=DEFAULT_CODES)
            domain.load_terms(db,term_table=term_table,mygo_database=opts[MYGO_DATABASE])
            domain.structures = structs
            domain.scale_structures()
            yield domain
  

    #_print("Loaded %i Domains" % count)                

def _structures(domain_key, fr=False, mcm=False, psi=False, query=None):
    """Load Structures for a domain key
    @return: Generator of (sf,confidence,domain_type)"""
    assert(fr | mcm | psi)
    f,m,p = (FRDomain, MCMDomain, PSIDomain)
    keep = {f:fr,m:mcm,p:psi}
    ginzu = {"fold_recognition":f, "msa":m, "pfam":m, "psiblast":p, "unassigned":m, "user_defined":m}
 
    db = _connection()
    cursor = db.cursor()
    query = query % domain_key
    _debug(domain_key, query)
    cursor.execute(query)
    for sccs,confidence,domain_type in cursor.fetchall():
            if not keep[ginzu[domain_type]] or sccs == None or sccs == "":
                continue
            _debug("\t", sccs, confidence, domain_type, keep[ginzu[domain_type]])
            # See if it's the type we want
            # Split the sccs, create domains and add the single structure
            for sf in sccs.split(","):
                sf = ".".join(sf.strip().split(".")[:3])
                _debug("\t","\t",sf,ginzu[domain_type])
            
                # I've mapped domain type strings to domain constructors
                # in the dict 'ginzu'.  This should be a factory!
                yield (sf,confidence,domain_type)
    cursor.close()



def sequences(sequence_keys=None, experiment_keys=None):
    """Return a list of sequences based on the selected options"""
    global _db
    #assert((sequence_keys != None) | (experiment_keys != None))
    db = _connection(db=opts[UPLOAD_DATABASE])
    cursor = db.cursor()
    table = opts[UPLOAD_TABLE]#"bayes"+opts[FUNCTION_SUFFIX]
    #kdrew: figure out a way to get parameterize domain_sccs
    if opts[SEQUENCE_TABLE]:
        query = "select distinct d.domain_sequence_key, d.parent_sequence_key from %s p join hpf.domain d on p.sequence_key=d.parent_sequence_key join %s s on d.domain_sequence_key=s.domain_sequence_key and s.sccs is not NULL left outer join %s b on s.domain_sequence_key=b.domain_sequence_key where (b.timestamp is NULL or b.timestamp < %s)" % (opts[SEQUENCE_TABLE],opts[DOMAIN_SCCS_TABLE], table,opts[OLDER])
    elif sequence_keys:
	print sequence_keys
        query = "select distinct s.domain_sequence_key, s.parent_sequence_key from %s s left outer join %s b on s.domain_sequence_key=b.domain_sequence_key where s.domain_sequence_key in (%s) and s.sccs is not NULL and (b.timestamp is NULL or b.timestamp < %s)" % (opts[DOMAIN_SCCS_TABLE],table,",".join(sequence_keys),opts[OLDER])
    elif experiment_keys:
        query = "select distinct d.domain_sequence_key, d.parent_sequence_key from hpf.protein p join hpf.domain d on p.sequence_key=d.parent_sequence_key join %s s on d.domain_sequence_key=s.domain_sequence_key and s.sccs is not NULL left outer join %s b on s.domain_sequence_key=b.domain_sequence_key where p.experiment_key in (%s) and (b.timestamp is NULL or b.timestamp < %s)" % (opts[DOMAIN_SCCS_TABLE], table,",".join([str(key) for key in experiment_keys]),opts[OLDER])
    else:
        query = "select distinct s.domain_sequence_key,s.parent_sequence_key from %s s left outer join %s b on s.domain_sequence_key=b.domain_sequence_key where s.sccs is not NULL and (b.timestamp is NULL or b.timestamp < %s)" % (opts[DOMAIN_SCCS_TABLE], table,opts[OLDER])
    _print(query)
    cursor.execute(query)
    #seqs += [(parent_key, domain_key) for id,sequence in cursor.fetchall()]
    seqs = [(pkey, dkey) for dkey, pkey in cursor.fetchall()]
    cursor.close();
    
    _print("Loaded %i sequences" % len(seqs))
    #_debug(seqs)
    _db.close()
    _db = None
    return seqs

def _connection(db=None):
    global _db
    if _db:
        return _db
    # Otherwise create a new connection
    if db:
        _db = MySQLdb.connect(db=db,read_default_file="~/.my.cnf")
    else:
        _db = MySQLdb.connect(read_default_file="~/.my.cnf")
        #_db = MySQLdb.connect(host="127.0.0.1", user="pwinters", passwd="bonneaulab", port=13307,db="hpf")
    return _db

def _datetime(time_str):
    return datetime.strptime(time_str, "%d %b %y")

def _mysqldate(dt):
    mysqldate = datetime.strftime(dt,"%Y-%m-%d %H:%M:%S")
    return mysqldate
    
def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "?duqe:s:mfpk:a:o:g:h:z:w:cbt:j:yx:n:", [DEBUG,UPLOAD,SILENT,EXPERIMENT,SEQUENCES,MCM,FR,PSI,SHELVE])
    global opts
    for o,a in _opts:
        if o in ('-?',"help"):
            _usage()
            sys.exit()
        if o == '-d':
            opts[DEBUG] = True
	#kdrew: was -x but conflicted with another option so changed to -u
        if o == '-u':
            opts[UPLOAD] = False
        if o == '-q':
            opts[SILENT] = True
        if o == '-e':
            opts[EXPERIMENT] = [int(e) for e in a.strip().split()]
        if o == '-s':
            assert(os.path.exists(a))
            opts[SEQUENCES] = a
        if o == '-m':
            opts[MCM] = True
        if o == '-f':
            opts[FR] = True
        if o == '-p':
            opts[PSI] = True
        if o == '-k':
            assert(os.path.exists(a))
            opts[SHELVE] = a
        if o == '-a':
            opts[PROCESSORS] = int(a)
        if o in ('-o',OLDER):
            opts[OLDER] = "DATE('"+_mysqldate(_datetime(a.strip()))+"')"
        if o in ('-g',FUNCTION_DATABASE):
            opts[FUNCTION_DATABASE] = a
        if o in ('-h',UPLOAD_DATABASE):
            opts[UPLOAD_DATABASE] = a
        if o in ('-z',FUNCTION_SUFFIX):
            opts[FUNCTION_SUFFIX] = "_"+a
        if o in ('-w',SEQUENCE_TABLE):
            opts[SEQUENCE_TABLE] = a
        if o in ('-n',DOMAIN_SCCS_TABLE):
            opts[DOMAIN_SCCS_TABLE] = a
        if o in ('-c',COMBINE):
            opts[COMBINE] = True
        if o in ('-b',BEST):
            opts[BEST] = True
        if o in ('-t',TERM_SUFFIX):
            opts[TERM_SUFFIX]=a
        if o in ('-j',UPLOAD_TABLE):
            opts[UPLOAD_TABLE]=a
        if o in ('-y',MCMDATA):
            opts[MCMDATA]=True
        if o in ('-x',STRUCTURE_QUERY):
            opts[STRUCTURE_QUERY] = a
	    print "query: "+a

    if opts[MCMDATA]:
        assert not (opts[PSI] or opts[FR]), "Using mcmdata table we can only process mcm domains"

    if not opts[TERM_SUFFIX]:
        opts[TERM_SUFFIX]=opts[FUNCTION_SUFFIX]
    
    # Table to Upload To
    if opts[UPLOAD_TABLE]:
        if not len(opts[UPLOAD_TABLE].split(".")) > 1:
            opts[UPLOAD_TABLE] = opts[UPLOAD_DATABASE]+"."+opts[UPLOAD_TABLE]
    else:
        opts[UPLOAD_TABLE] = opts[UPLOAD_DATABASE]+".bayes"+opts[FUNCTION_SUFFIX]
    
    assert(opts[MCM] | opts[PSI] | opts[FR])
    return opts, args

def _debug(*args):
    if opts[DEBUG]:
        _print(*args)

def _print(*args):
    if not opts[SILENT]:
        print " ".join([str(a) for a in args])

def _usage():
    print "$ function_prediction.py [-options]"
    print "Perform Function Prediction for a given set of sequences or experiment."
    print ""
    print "Options"
    print "-d debug mode"
    print "-q silent"
    print ""
    print "-e experiment number to grab sequences for."
    print "-s file to load sequence keys from."
    print "-w sequence table"
    print ""
    print "-m use mcm results"
    print "-f use fold_recognition results"
    print "-p use psiblast results"
    print ""
    print "-c combine superfamilies"
    print "-b choose best superfamily score for domain"
    print "-y use the mcmdata table for structure scores"
    print "-x Structure Query, should return a tuple (sccs,confidence,domain_type)"
    print ""
    print "-h database to upload to"
    print "-z function table suffix to upload to, ie '_goLite_062009_3'"
    print "-g database that holds function tables"
    print "-j table to upload results"
    print "-t suffix for the iea term table"
    print ""
    print "-o update domains with records older than this date"
    print "     Format: '30 Jun 01' from time.strptime(opt, '%d %b %y')"
    print "-k load/store dictionary shelves from dir argument. expects 'LogRatios_*.shelve', 'MutualInformation_*.shelve'"
    
if __name__=="__main__":
    _do(sys.argv[1:])


