#!/usr/bin/env python

import numpy
import sys
import getopt
from hpf.function.metric import Metric
from hpf.processing import MultiProcessor
import MySQLdb
import shelve

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
SHELVE = "shelve"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, 
        SHELVE:"/Users/kdrew/data/hpf/functionTables/probability_goLite_062009.shelve"}

__hpf = None

_prob = None
_mf_acc = None

def calc_corr(sf1,sf2):
    score1 = []
    score2 = []
    weights = []
    union = 0
    intersection = 0
    for mf in _mf_acc:
        # P(mf|sf)
        a=_prob.get_metric(sf1,mf) 
        b=_prob.get_metric(sf2,mf)
        #print sf1,sf2,mf,a,b
        if a==0 and b==0:
            continue
        if a>0 and b>0:
            intersection += 1
        union+=1
        w=_prob.get_metric(mf)
        # Append the values to the appropriate list
        for v,l in [(a,score1),(b,score2),(w,weights)]:
            l.append(v)

    assert len(score1)==len(score2)
    assert len(score1)==len(weights)
    assert len(weights)==len(score2)
    if len(score1) == 0 or intersection == 0:
        return

    score1 = numpy.array(score1)
    score2 = numpy.array(score2)
    weights = numpy.array(weights)
    
    corr_coeff = corr(score1,score2,weights)
    if corr_coeff == numpy.nan:
        return
    return (sf1, sf2, corr_coeff, intersection, union)
    
def __upload(result):
    if result==None:
        return
    sf1, sf2, corr_coeff, intersection, union = result
    
    db = _hpf(db="pdb")
    cursor = db.cursor()
    query = """
        insert into correlation 
        (sf1, sf2, corr_coeff, intersection, union_mf)
        VALUES (%s,%s,%s,%s,%s)
        on duplicate key update 
        sf1=values(sf1),
        sf2=values(sf2),
        corr_coeff=values(corr_coeff),
        intersection=values(intersection),
        union_mf=values(union_mf),
        timestamp=NOW()
        """
    cursor.executemany(query,[(sf1,sf2,corr_coeff,intersection,union)])
    cursor.close()
    #print result
    
def __calc_corr(task):
    sf1, sf2 = task
    return calc_corr(sf1,sf2)

def main():
    # Get all superfamilies and molecular functions
    global _mf_acc
    with MySQLdb.connect(db="functionTables",passwd="patrick_nyu") as cursor:
        # Only use superfamilies with something in the probability table
        query = """
            select distinct acc 
            from functionTables.probability_goLite_062009 
            where acc like '%.%'
            """
#        query = """select distinct substring_index(sccs,'.',3) 
#            from pdb.astral95_1_75 a 
#            join functionTables.probability_golite_062009 p
#            on substring_index(a.sccs,'.',3)=p.acc"""
        print query
        cursor.execute(query)
        sccs = [t[0] for t in cursor.fetchall()]
        print len(sccs)," superfamilies"
        # Only use molecular functions
        query = """
            select distinct p.acc 
            from functionTables.probability_goLite_062009 p 
            join mygoLite_062009.term t 
            on p.acc=t.acc and t.term_type='molecular_function' and t.acc!='GO:0003674'
            where p.acc2 is NULL
            order by p.metric asc
            """
        print query
        cursor.execute(query)
        _mf_acc = [t[0] for t in cursor.fetchall()]
        print len(_mf_acc)," molecular functions"

    tasks=[]
    for i,sf1 in enumerate(sccs):
        for sf2 in sccs[i+1:]:
            tasks.append((sf1,sf2))
    print len(tasks)," pairwise superfamilies"
            
    global _prob
    print "Opening shelve"
    dict = shelve.open(opts[SHELVE])
    _prob = Metric(dict=dict,default=0)
    pool = MultiProcessor(8, modulus=100, raise_errors=True)
    for r in pool.run(__calc_corr, tasks, __upload):
        pass

def m(x,w):
    return sum(x*w)/sum(w)

def cov(x,y,w):
    return sum( w*(x-m(x,w))*(y-m(y,w)) )/sum(w)

def corr(x,y,w):
    return cov(x,y,w)/numpy.sqrt(cov(x,x,w)*cov(y,y,w))

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
