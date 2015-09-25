#!/usr/bin/env python
import sys
import MySQLdb
import getopt

PROTEIN = "protein"
TAX = "tax"
EXPERIMENT = "experiment"
LENGTH = "length"
TUNNEL = "tunnel"
opts = {PROTEIN: True, TAX: True, EXPERIMENT: None, LENGTH:False, TUNNEL:False}


def proteins(e_ids):
    query = "select s.sequence_key, s.sequence, length(s.sequence), group_concat(gi) from (select distinct p.sequence_key as sequence_key, s.sequence as sequence from hpf.experiment e join hpf.protein p join hpf.sequence s on e.id=p.experiment_key and p.sequence_key=s.id where e.id in (%s)) s left outer join hpf.sequenceAc a on s.sequence_key=a.sequence_key group by s.sequence_key" % e_ids
    return query
        
def domains(e_ids):
    query = "select s.sequence_key, s.sequence, length(s.sequence), group_concat(gi) from (select distinct d.domain_sequence_key as sequence_key, s.sequence as sequence from hpf.experiment e join hpf.protein p join hpf.domain d join hpf.sequence s on e.id=p.experiment_key and p.sequence_key=d.parent_sequence_key and d.domain_sequence_key=s.id where e.id in (%s)) s left outer join ddbCommon.ddb_sequenceAc a on s.sequence_key=a.sequence_key group by s.sequence_key" % e_ids
    return query

def main():
    if opts[PROTEIN]:
        func = proteins
    else:
        func = domains

    db = _connect()
    cursor = db.cursor()
    cursor.execute(_e_name())
    for id,ename,nname,ncat in cursor.fetchall():
        if nname:
            if ncat != 'scientific name':
                continue
            name = nname
        else:
            name = ename
        name = name.replace(" ","_")
        name = name.replace(".","_")
        name = name.replace("/","_")
        name = name.replace("\\","_")
        
        type = (opts[PROTEIN] and ["p"] or ["d"])[0]
        out = open("%s.%i.%s.fa" % (name,id,type), "w")
        # Now get entries
        cursor.execute(func(str(id)))
        print id,name
        for key, sequence, length, gis in cursor.fetchall():
            if gis:
                entry = ">%s|%s\n%s" % (key, gis, sequence.upper())
            else:
                entry = ">%s\n%s" % (key, sequence.upper())
            print >>out, entry
    out.close()
    cursor.close(); db.close()

def _e_name():
    if opts[EXPERIMENT]:
        where = " where e.id=%s"%opts[EXPERIMENT]
    else:
        where = ""
    
    if opts[TAX]:
        return "select e.id,e.name,n.name,n.category from hpf.experiment e left outer join ddbCommon.nrTaxName n on e.taxonomy_id=n.taxonomy_id"+where
    else:
        return "select e.id,e.name,NULL,NULL from experiment e"+where

def _connect():
    if opts[TUNNEL]:
        return MySQLdb.connect(host="127.0.0.1", db="hpf", user="pwinters", passwd="bonneaulab", port=3307)
    else:
        return MySQLdb.connect(host="localhost", db="hpf", user="patrick", passwd="patrick_nyu", port=3306)

def _do(argv):
    try:
        opts, args = _options(argv)
    except Exception, e:
        print e
        _usage()
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "dne:l")
    global opts
    for o,a in _opts:
        if o == '-d':
            opts[PROTEIN] = False
        if o == '-n':
            opts[TAX] = False
        if o == '-e':
            opts[EXPERIMENT] = int(a)
        if o == '-l':
            opts[LENGTH] = True
        if o == '-t':
            opts[TUNNEL] = True
    return opts, args
    
def _usage():
    print "$ template.py [-options] args"
    print "Dump sequences into fasta files based on experiment."
    print ""
    print "Options"
    print "-d print domain sequences as opposed to protein sequences"
    print "-n name experiment files based on experiment name rather than scientific taxonomy name"
    print "-e integer experiment id to dump out"
    print "-l print length of sequence with identifier"
    sys.exit(2)

if __name__=="__main__":
    _do(sys.argv[1:])

