#!/usr/bin/env python
import sys
import MySQLdb
from script.util import conditional
from collections import defaultdict


print "Parsed records"



#!/usr/bin/env python

import sys
import getopt

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "upload"
XML = "xml"
TABULAR = "tabular"
ALIGN = "alignment"
EXPECT = "expect"
DATABASE = "database"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, XML:True, TABULAR:False, EXPECT:1e-5, ALIGN: 0.75 ,DATABASE:"since_solved"}

def text_records(input):
    from Bio.Blast import NCBIStandalone
    blast_parser = NCBIStandalone.BlastParser()
    blast_iterator = NCBIStandalone.Iterator(input, blast_parser)
    for blast_record in blast_iterator:
        yield blast_record
        
def xml_records(input):
    from Bio.Blast import NCBIXML, NCBIStandalone
    _debug("Parsing",input)
    for blast_record in NCBIXML.parse(input):
        yield blast_record

def _upload(cursor, query, alignment_id, expect,align_perc, experiment_name):
    """Expects a records iterable or generator"""
    db_query = """
        insert ignore into blast 
        (query, subject, expect, align_perc, experiment) 
        values ('%s','%s',%e,%f,'%s')
        """ % (query, alignment_id, expect,align_perc,experiment_name)
    cursor.execute(db_query)

def _connect(db=None, host="localhost"):
    _db = MySQLdb.connect(db=db,read_default_file="~/.my.cnf")
    return _db



def tabular(input, upload=False, expect_cutoff = 1e-3, align_cutoff = 0.75, experiment_name=None):
    """Process tabular blast ouput"""
    def _threshold(align_perc, expect):
        """Defines a threshold for an HSP to pass"""
        return float(align_perc) >= align_cutoff and float(expect) <= expect_cutoff  
    
    scores = defaultdict(list)
    
    for line in input:
        # Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
        if line.startswith("#"):
            continue
        query_id, subject_id, percent_ident, align_len, mismatches, gaps, q_start, q_end, s_start, s_end, expect, score = line.strip().split("\t")
        tuple = (expect, subject_id, align_len, experiment_name)
        scores[query_id].append( tuple )
        
    ids = ",".join([sequence_key for sequence_key in scores])
    print len(ids)
    db = _connect(db=opts[DATABASE])
    cursor = db.cursor()
    query = "select distinct s.id, length(s.sequence) from sequence s where s.id in (%s)" % ids
    _debug(query)
    cursor.execute(query)
    records = []
    sequence_ids = set()
    for id,length in cursor.fetchall():
        for expect, subject_id, align_len, experiment_name in scores[str(id)]:
            align_perc = float(align_len)/float(length)
            if _threshold(align_perc, expect):
                tuple = (id,subject_id,expect,align_perc,experiment_name)
                records.append(tuple)
                #_print(tuple)
                sequence_ids.add(id)
    #print ", ".join([str(id) for id in sequence_ids])
    print len(set(sequence_ids))
    

def records(iterator, upload=False, expect_cutoff = 1e-3, align_cutoff = 0.75, experiment_name=None):
    """Iterate over records, perform thresholding, and upload if necessary"""   
    def _threshold(align_perc, expect):
        """Defines a threshold for an HSP to pass"""
        return float(align_perc) >= align_cutoff and float(expect) <= expect_cutoff  
    
    scores = defaultdict(list)
    done=0
    if upload:
        db = _connect(db=opts[DATABASE])
        cursor = db.cursor() 
    
    _debug("Iterating")
    for blast_record in iterator:
        _debug(blast_record)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
               
                length_of_longer_sequence = max(blast_record.query_letters, alignment.length)
                length_of_alignment = max(len(hsp.query),len(hsp.match))
                align_perc = float(length_of_alignment)/float(length_of_longer_sequence)
                
                if _threshold(align_perc, hsp.expect):
                    scores[blast_record.query].append((alignment,hsp))
                    alignment_id = alignment.title[1:].split()[0]
                    if upload:
                        _upload(cursor, blast_record.query, alignment_id,hsp.expect,align_perc,experiment_name)
        
        scores[blast_record.query].sort(cmp=lambda x,y: cmp(x[1].expect,y[1].expect))
        #records.append(blast_record)
        done+=1
        #if done % 1000 == 0:
        _print("Finished parsing",blast_record.query,done,"done")
    return scores

def main(experiment_name, input=None):
    if input == None:
        input = sys.stdin
    elif isinstance(input, str):
        input = open(input)
    if opts[TABULAR]:
        tabular(input, expect_cutoff=opts[EXPECT], align_cutoff=opts[ALIGN], experiment_name=experiment_name)
    else:
        records(conditional(opts[XML], xml_records, text_records)(input), upload=opts[UPLOAD], expect_cutoff=opts[EXPECT], align_cutoff=opts[ALIGN], experiment_name=experiment_name)
    
def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "dxqe:a:th:")
    assert(len(args)>=1)
    global opts
    for o,a in _opts:
        if o == '-d':
            opts[DEBUG] = True
        if o == '-u':
            opts[UPLOAD] = True
        if o == '-x':
            opts[XML] = True
        if o == '-t':
            opts[TABULAR] = True
        if o == '-q':
            opts[SILENT] = True
        if o == '-e':
            opts[EXPECT] = float(a)
        if o == '-a':
            opts[ALIGN] = float(a)
        if o == '-h':
            opts[DATABASE] = a
    return opts, args

def _debug(*args):
    if opts[DEBUG]:
        _print(*args)

def _print(*args):
    if not opts[SILENT]:
        print " ".join([str(a) for a in args])
    
def _usage():
    print "$ blast_load.py [-options] experiment_name [input_file]"
    print "Parse and/or upload Blast results from input_file or stdin"
    print ""
    print "Options"
    print "-d Debug"
    print "-q Quiet"
    print "-x XML format"
    print "-u Upload results"
    print "-e Expectation Cutoff"
    print "-a Alignment Percentage of Larger Sequence Cutoff"
    print "-h Database to upload to."


if __name__=="__main__":
    _do(sys.argv[1:])


