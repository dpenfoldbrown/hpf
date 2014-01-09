#!/usr/bin/env python

import MySQLdb
import sys
import re

input = sys.stdin
db_name = sys.argv[1]

db = MySQLdb.connect(user="patrick", passwd="patrick_nyu", db="since_solved")
cursor = db.cursor()

cluster=0
cluster_entry_num=0
for line in input:
    try:
        if line.startswith(">Cluster"):
            cluster = int(line.strip().split()[1])
        else:
            line = line.strip()
            end = line.find("|")
            if end==-1:
                end=line.find("...")
            if end==-1:
                end=len(line)
            identifier = line[line.index(">")+1:end]
            cluster_entry_num = int(line.split()[0])
            match = re.compile("[0-9]+aa").search(line)
            aa_count = int(line[match.start():match.end()][:-2])    
            # execute SQL statement
            cmd = "INSERT INTO cdhit_clstr (cluster_id, cluster_entry_num, aa_count, identifier, db, foreign_key) values (%i, %i, %i, '%s', '%s', %s)" % (cluster, cluster_entry_num, aa_count, identifier, db_name, identifier)
            #print cmd
            cursor.execute(cmd)
    except:
        print line
        raise
