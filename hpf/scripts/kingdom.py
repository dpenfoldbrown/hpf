#!/usr/bin/env python
import MySQLdb
from collections import defaultdict

def _connection():
    return MySQLdb.connect(host="127.0.0.1", user="pwinters", passwd="bonneaulab", db="hpf", port=13307)    

kingdoms = {2:"Bacteria",2157:"Archaea",2759:"Eukaryota"}

def tax_recurse(tax_id, cursor):
    if tax_id in kingdoms.keys():
        return tax_id
    else:
        cursor.execute("select parent_taxonomy_id from ddbCommon.nrTaxNode where taxonomy_id=%i" % tax_id)
        parent = cursor.fetchone()[0]
        #print "\t",tax_id,parent
        return tax_recurse(parent,cursor)


con = _connection()
query = "select id, taxonomy_id from hpf.experiment"
cursor = con.cursor()
cursor.execute(query)

all = defaultdict(lambda: list())
for eid,taxid in cursor.fetchall():
    try:
        kingdom_key = tax_recurse(taxid,cursor)
    except:
        continue
    print eid, kingdoms[kingdom_key]
    all[kingdom_key].append(eid)
    
for kingdom_key in all:
    print kingdoms[kingdom_key]
    print ", ".join([str(e) for e in all[kingdom_key]])
 
