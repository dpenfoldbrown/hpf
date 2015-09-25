#!/usr/bin/env python

import sys, os
import getopt
import MySQLdb
from hpf.runtime import runtime, Flag, IntegerOption
from hpf.utilities import consume
from hpf.processing import processor,SYNCHRONOUS, PROCESSORS
from Bio.Entrez import efetch,esearch, read as eread

COUNT = "Count"
IDLIST = "IdList"
MAX_RETURN = 200


def tasks():
    query = "select distinct gi from domain_sccs d join sequenceAc a on d.parent_sequence_key=a.sequence_key where ac is not NULL and d.sccs is not NULL and d.domain_type !='psiblast'"
    with MySQLdb.connect(host="127.0.0.1",passwd="patrick_nyu",db="hpf") as cursor:
        runtime().debug(query)
        cursor.execute(query)
        runtime().debug("Fetching")
        all = cursor.fetchall()
        return [gi[0] for gi in all]

def _publications(gi):
    return publications(gi)

def pubmed(gi,ids,query):
    """
    Get the pubmed articles listed by *ids
    """
    _ids=",".join(ids)
    for id in ids:
        handle = efetch(db="pubmed",id=id,retmode='xml',rettype='xml',retmax=MAX_RETURN)
        try:
            #print handle.read()
            results = eread(handle)
            for citation in results:
                #runtime().debug(citation.keys())
                citation = citation['MedlineCitation']
                pmid = citation['PMID']
                article = citation['Article']
                title = article['ArticleTitle']
                journal = article['Journal']['Title']
                try:
                    date = citation['DateCompleted'] if citation.has_key('DateCompleted') else citation['DateCreated']
                    year = date['Year']
                    month = date['Month']
                    day = date['Day']
                    datetime = "%s-%s-%s" % (year,month,day)
                except:
                    datetime = '0000-00-00'
                
                runtime().debug("Parsed pmid:%s" % id)
                yield Citation(gi, pmid, title, journal, datetime, query)
        except:
            runtime().debug("Failure fetching pmid:%s" % id)
            continue
        finally:
            handle.close()

class Citation(object):
    def __init__(self,gi,pmid, title,journal,date,query):
        self.gi=gi
        self.pmid = pmid
        self.title=title
        self.journal=journal
        self.date=date
        self.query=query
        
    def __repr__(self):
        #return str(self.gi)+str(self.pmid)
        return "%i %s %s %s" % (self.gi,
                                   self.date, 
                                   self.title, 
                                   self.pmid)
    
    def _tuple(self):
        return (self.gi,self.query,self.pmid,self.title,self.journal,self.date)

def upload(citations):
    if citations == None:
        return
    query = """INSERT INTO hpf.publications (gi, query, pmid, title, journal, date)
        VALUES (%s, %s, %s, %s, %s, %s) 
        ON DUPLICATE KEY UPDATE 
            gi=VALUES(gi),
            query=VALUES(query),
            pmid=VALUES(pmid),
            title=VALUES(title),
            journal=VALUES(journal),
            date=VALUES(date)"""
    with MySQLdb.connect(host="127.0.0.1",db="hpf",passwd="patrick_nyu") as cursor:
        for cit in citations:
            try:
                cursor.executemany(query, [cit._tuple()])
                runtime().debug("\t",str(cit))
            except:
                continue


def names(gi):
    columns = ["protgi","gi"]
    tables = [("accession",0),
              #("bind",0),
              ("codedby",0),
              ("genbank",0),
              #("geneid",0),
              ("genename",0),
              #("kegg",0),
              ("locustag",0),
              ("pdb",0),
              ("pfam",0),
              ("unigene",0),
              ("uniprotkb",0)]
    query = "select * from %s where %s=%i"
    all = set()
    with MySQLdb.connect(db="synonyms3",host="err.bio.nyu.edu",user="patrick",passwd="patrick_nyu") as cursor:
        for name, col in tables:
            table = "refseq_%s" % name
            column = columns[col]
            q = query % (table,column,gi)
            runtime().debug(q)
            cursor.execute(q)
            rows = cursor.fetchall()
            for row in rows:
                if row[1]!=None:
                    all.add(str(row[1]))
    return all
            

def publications(gi):
    """
    Grab the list of publications for a given protein
    """
    gi_names = list(names(gi))
    if len(gi_names) == 0:
        return
    #runtime().debug("NAMES",gi_names)
    term = " OR ".join(['"%s"' %n for n in gi_names])
    runtime().debug(term)
    handle = esearch(db="pubmed", 
            term=term,
            retmax=MAX_RETURN,
            mindate="2000")
    try:
        results = eread(handle)
    finally:
        handle.close()
    results[COUNT] = int(results[COUNT])
    runtime().debug("Found",results[COUNT])
    if results[COUNT]>0:
        runtime().debug(gi,term,"\t",results[COUNT], "\n\t",results[IDLIST])
        return list(pubmed(gi,set(results[IDLIST]),term))
        

def main(*args):
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS),processors=runtime().opt(PROCESSORS))
    runtime().debug("Using processor",pool)
    pool.make_tasks(tasks)
    consume(pool.run(_publications,result=upload))
    
def _do(argv):
    r = runtime()
    r.description("""
    publications.py [-options] args
    Gather the number and type of publications for a gi.
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing"))
    r.add_option(IntegerOption(PROCESSORS, 'a', description="Number of processors for running on"))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
