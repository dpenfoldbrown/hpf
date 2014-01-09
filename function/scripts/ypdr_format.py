#!/usr/bin/env python
#Formats the predictions for the YPDR
#Adjusts scores according to error model
#IF NO FORMATTING JUST DO A QUERY

raise Exception("Don't run this script.  It's no longer necessary.  You can upload results directly from mysql because we aren't formatting with error models.")

import os
import sys
import getopt
import MySQLdb
from hpf.processing import MultiProcessor

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
opts = {DEBUG: False, SILENT: False, UPLOAD: True}

BAYES_TABLE = "functionPredictions.bayes_golite_062009_3"

_db = None
_cursor = None

eukaryotes = set([804, 822, 823, 825, 826, 827, 886, 915, 924, 900, 888, 917, 4])
prokaryotes = set([29, 30, 31, 34, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 889, 920])

# CARL bayes_llr_results
#+---------------------+------------------------+------+-----+-------------------+-------+
#| Field               | Type                   | Null | Key | Default           | Extra |
#+---------------------+------------------------+------+-----+-------------------+-------+
#| parent_sequence_key | int(11)                | YES  | MUL | NULL              |       |
#| domain_sequence_key | int(11)                | NO   | PRI | 0                 |       |
#| mf_acc              | varchar(20)            | NO   | PRI |                   |       |
#| name                | varchar(255)           | NO   |     |                   |       |
#| p_llr               | double                 | YES  |     | 0                 |       |
#| l_llr               | double                 | YES  |     | 0                 |       |
#| s_llr               | double                 | YES  |     | 0                 |       |
#| pl_llr              | double                 | YES  |     | 0                 |       |
#| ps_llr              | double                 | YES  |     | 0                 |       |
#| ls_llr              | double                 | YES  |     | 0                 |       |
#| pls_llr             | double                 | YES  |     | 0                 |       |
#| base_llr            | double                 | YES  |     | NULL              |       |
#| type                | enum('mcm','fr','psi') | YES  | MUL | NULL              |       |
#| timestamp           | timestamp              | NO   | MUL | CURRENT_TIMESTAMP |       |
#+---------------------+------------------------+------+-----+-------------------+-------+

def eukaryote(experiment_key):
    return experiment_key in eukaryotes
    
def prokaryote(experiment_key):
    return experiment_key in prokaryotes

def mcm(experiment_key):
    assert(false)
    if prokaryote(experiment_key):
        return [5.25,3.5,2.7]
    elif eukaryote(experiment_key):
        return [6.7,8.5,4.0]
    
def fr(experiment_key):
    assert(false)
    if prokaryote(experiment_key):
        return [3.8,4.5,1.9]
    elif eukaryote(experiment_key):
        return [25.0,16.0,12.0]
    
def psi(experiment_key):
    assert(false)
    if prokaryote(experiment_key):
        return [1.3,3.0,5.0]
    elif eukaryote(experiment_key):
        return [25.0,10.0,5.0]
    
def bin_llr(base_llr):
    if base_llr >= -4:
        return 0
    elif base_llr >= -8:
        return 1
    else:
        return 2



# BDDB GO
#
#                 id: 2211763
#       sequence_key: 1658245
#domain_sequence_key: 1658245
#                acc: GO:0004930
#               name: G-protein coupled receptor activity
#          term_type: molecular_function
#      evidence_code: KD
#        xref_dbname: 
#           xref_key: 
#        probability: 0
#                llr: 2.27587478562
#              level: 0
#             source: kd826_psi_bayes8_s_1
#        insert_date: NULL
#          timestamp: 2008-07-16 11:46:58
def format(row):
    domain_type = row["type"]
    #experiment_key = row["experiment_key"]
    #models = {"mcm":mcm,"fr":fr,"psi":psi}
    #model = models[domain_type](experiment_key)
    #if model==None:
    #    raise Exception("Unknown kingdom",experiment_key)
    #bin = bin_llr(row["base_llr"])
    #threshold = model[bin]
    pls_llr = row["pls_llr"]#-threshold
    source = domain_type+"_bayes_pls_golite062009"
    d = {"sequence_key":row["parent_sequence_key"],
         "domain_sequence_key":row["domain_sequence_key"],
         "acc":row["acc"],
         "name":row["name"],
         "term_type":"molecular_function",
         "evidence_code":"KD",
         "probability":0,
         "llr":pls_llr,
         "level":0,
         "source":source,
         "insert_date":row["timestamp"],
         "xref_key":"",
         "xref_dbname":""}
    return d
    
def upload(result):
    # Upload everything and we can threshold later to limit what gets in Lars' db.
    if opts[UPLOAD]:# and result["llr"] > -10.0:
        #print result
        cursor = _bddb()
        cursor.executemany("""Insert into hpf.go 
            (sequence_key, 
             domain_sequence_key, 
             acc, 
             name, 
             term_type,
             evidence_code, 
             probability,
             llr,level,
             source,
             insert_date,
             xref_dbname,
             xref_key) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
             ON DUPLICATE KEY UPDATE 
                probability=VALUES(probability),
                llr=VALUES(llr),
                source=VALUES(source),
                insert_date=VALUES(insert_date)
                """
            ,[(result["sequence_key"],result["domain_sequence_key"],result["acc"],
               result["name"],result["term_type"],result["evidence_code"],
               result["probability"],result["llr"],result["level"],
               result["source"],result["insert_date"],
               result["xref_dbname"],result["xref_key"])])

def main():
    func_db = _func()
    func_cursor = func_db.cursor(MySQLdb.cursors.DictCursor)
    query = """
        select b.parent_sequence_key,b.domain_sequence_key,
        b.mf_acc as acc, b.name, b.pls_llr,b.base_llr,b.type,b.timestamp 
        from %s b 
        where pls_llr > 0
        """ % BAYES_TABLE
    print query
    func_cursor.execute()
    tasks = func_cursor.fetchall()
    func_cursor.close()
    func_db.close()
    pool = MultiProcessor(processors=8, modulus=100, raise_errors=False)
    pool.run(format, tasks, upload)

def _func():
    return MySQLdb.connect(db="hpf",read_default_file="~/.my.cnf")
    
def _bddb():
    global _db, _cursor
    if not _db:
        _db = MySQLdb.connect(host="127.0.0.1",port=13307,user="pwinters",passwd="bonneaulab",db="hpf")
    if not _cursor:
        _cursor = _db.cursor()
    return _cursor

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
    print "$ ypdr_format.py [-options] args"
    print "Adjust, format, and upload function predictions into the bddb GO table."
    print ""
    print "Options"


if __name__=="__main__":
    _do(sys.argv[1:])
