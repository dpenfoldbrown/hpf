#!/usr/bin/env python
import MySQLdb
import sys

experiments = sys.argv[1:]
db = MySQLdb.connect(host="127.0.0.1", user="pwinters", passwd="bonneaulab", db="hpf", port=13307)
cursor = db.cursor()

query = "select distinct s.id,s.sequence from bddb.experiment e join bddb.protein p join ddbCommon.domain d join ddbCommon.ddb_sequence s join bddb.mcmData m on e.id=p.experiment_key and p.sequence_key=d.parent_sequence_key and d.domain_sequence_key=s.id and s.id=m.sequence_key"
if len(experiments) > 0:
    in_statement = " where e.id in (%s)" % ", ".join(experiments)
    query += "%s" % in_statement
#print query
cursor.execute(query)
for hpf_id,sequence in cursor.fetchall():
    print ">%s\n%s" % (hpf_id, sequence)
