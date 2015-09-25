#!/usr/bin/env python
import MySQLdb
import sys
import traceback

#kdrew: the sql command to create the table used in this script
#create table astral (id INT NOT NULL AUTO_INCREMENT PRIMARY KEY, INDEX(id), astral_id VARCHAR(32), pdb_id VARCHAR(32), sccs VARCHAR(32), chain_id VARCHAR(32), sequence longtext, timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP);


db = MySQLdb.connect(host="localhost", user="patrick", passwd="patrick_nyu", db="pdb")
cursor = db.cursor()

def fasta():
    cursor.execute("SELECT astral_id, pdb_id, sequence from astral")
    for record in cursor.fetchall():
        astral_id, pdb_id, sequence = record
        print ">%s|%s\n%s" % (astral_id, pdb_id, sequence)
    

def insert():
    astral_table = "astral95_1_75"
    astral_id = None
    pdb_id = None
    sequence = ""
    for line in sys.stdin:
        line = line.strip()
        if line.startswith(">") or len(line)==0:
            if astral_id:
                cursor.execute("INSERT INTO %s (astral_id, pdb_id, sccs, chain_id, sequence) VALUES ('%s','%s','%s','%s','%s')" % (astral_table, astral_id, pdb_id, sccs, chain_id, sequence))
            if len(line)==0:
                continue
            sequence = ""
            #example
            #>d1uvxa_ a.1.1.1 (A:) Protozoan/bacterial hemoglobin {Green alga (Chlamydomonas eugametos) [TaxId: 3053]}
            split = line.split()
            astral_id = split[0][1:]
            pdb_id = astral_id[1:5]
            sccs = split[1]
            chain_id = split[2][1:-1]
            comment = MySQLdb.escape_string(" ".join(split[3:]))
        else:
            sequence += line

if __name__=="__main__":
    insert()
