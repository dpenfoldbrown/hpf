#!/usr/bin/env python
#>gi|15230298| Arabidopsis thaliana [YPDR-SEQID 312487]
#MTASGGGSTAATGRMPTWKERENNKKRERRRRAIAAKIFTGLRSQGNYKLPKHCDNNEVLKALCLEAGWIVHEDGTTYRKGSRPTETTVPCSSIQLSPQSSAFQSPIPSYQASPSSSSYPSPTRFDPNQSSTYLIPYLQNLASSGNLAPLRISNSAPVTPPISSPRRSNPRLPRWQSSNFPVSAPSSPTRRLHHYTSIPECDESDVSTVDSCRWGNFQSVNVSQTCPPSPTFNLVGKSVSSVGVDVSVKPWEGEKIHDVGIDDLELTLGHNTKGRG

import sys, os
import getopt
import MySQLdb

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
GI = "gi"

opts = {DEBUG: False, SILENT: False, UPLOAD: True, GI: False}
_db=None

def main():
    connection = _connection()
    cursor = connection.cursor()
    query = """
    Select distinct p.sequence_key, a.gi, s.sequence, n.name, e.name
    from hpf.experiment e join hpf.protein p on e.id=p.experiment_key
    join hpf.sequence s on p.sequence_key=s.id
    left outer join hpf.nrTaxName n on e.taxonomy_id=n.taxonomy_id
    left outer join hpf.sequenceAc a on s.id=a.sequence_key
    where n.category='scientific name' or n.category is NULL
    group by p.sequence_key
    """ 
    _debug(query)
    cursor.execute(query)
    for seq_key,gi,seq,tax_name,exp_name in cursor.fetchall():
        if tax_name:
            name = tax_name
        else:
            name = exp_name
        if gi and opts[GI]:
            print ">gi|%s| %s [BDDB-SEQID %s] \n%s" % (gi,name,seq_key,seq)
        else:
            print ">bddb|%s| %s [BDDB-SEQID %s] \n%s" % (seq_key,name,seq_key,seq)

def _connection(db=None):
    global _db
    if not _db:
        _db = MySQLdb.connect(read_default_group="~/.my.cnf")
    return _db

def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "dxqg", [DEBUG,SILENT,UPLOAD,GI])
    global opts
    for o,a in _opts:
        if o in ('-d',DEBUG):
            opts[DEBUG] = True
        if o in ('-g',GI):
            opts[GI] = True
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
