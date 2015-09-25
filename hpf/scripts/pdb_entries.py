#!/usr/bin/env python
#ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx
import sys
import MySQLdb

db = MySQLdb.connect(host="localhost", user="patrick", passwd="patrick_nyu", db="pdb")
cursor = db.cursor()

#Skip the first two headers
sys.stdin.readline()
sys.stdin.readline()

for line in sys.stdin:
    # Declare
    IDCODE, HEADER, ASCESSION_DATE, COMPOUND, SOURCE, AUTHOR_LIST, RESOLUTION, EXPERIMENT_TYPE = (None,None,None,None,None,None,None,None)
    parts = [IDCODE, HEADER, ASCESSION_DATE, COMPOUND, SOURCE, AUTHOR_LIST, RESOLUTION, EXPERIMENT_TYPE]
    split = line.strip().split("\t")
    for i in range(len(parts)):
        if i==2:
            month, day, year = split[i].split("/")
            parts[i] = "-".join([year,month,day])
        else:
            parts[i] = MySQLdb.escape_string(split[i])
        
    
    query = "insert into pdb_entries (idcode, header, ascession_date, compound, source, author_list, resolution, experiment_type) values ('%s', '%s', DATE('%s'), '%s', '%s', '%s', '%s', '%s');" % tuple(parts)
    cursor.execute(query)
