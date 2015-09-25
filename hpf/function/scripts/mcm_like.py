#!/usr/bin/env python
"""
Select a subset of psiblast domains that have superfamilies with a similar
distribution to MCM predictions.
"""

import os
import sys
import getopt
import MySQLdb

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
opts = {DEBUG: False, SILENT: False, UPLOAD: True}

#Selects percentage of each superfamily.  The 5227 is number of domains with MF.
query="""
select substring_index(sccs,'.',3), count(distinct d.domain_sequence_key)/5227 
from domain_sccs d join hddb_IEA.hddb_iea_golite_062009 g 
on d.parent_sequence_key=g.sequence_key and g.term_type='molecular_function' 
and g.acc!='GO:0003674' 
where domain_type in ('msa','pfam','unassigned') 
group by substring_index(sccs,'.',3) 
order by count(distinct d.domain_sequence_key) desc;
"""

def main():
    with _connection() as cursor:
        cursor.execute(query)
        all = cursor.fetchall()
        first = True
        select_query = """
select distinct d.domain_sequence_key, sccs, domain_type from domain_sccs d 
join hddb_IEA.hddb_iea_golite_062009 g 
on d.parent_sequence_key=g.sequence_key and g.term_type='molecular_function' 
and g.acc!='GO:0003674'
where sccs like '%s' and domain_type='psiblast' order by rand() limit %i""" 
        for sf,perc in all:
            number = int(5000*perc)
            if not number>0:
                continue
            #print sf,number
            prefix = "create table hpf_test.psi_mcm_like " if first else "insert into hpf_test.psi_mcm_like "
            q = prefix+select_query%("%"+sf+"%",number)
            print sf,number,q
            cursor.execute(q)
            first = False
            
    
def _connection():
    return MySQLdb.connect("localhost",user="patrick",passwd="patrick_nyu",db="hpf")

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
