#!/usr/bin/env python
import MySQLdb
import sys
import traceback

if __name__=="__main__":
    db = MySQLdb.connect(host="mcpeepants.bio.nyu.edu", user="patrick", passwd="patrick_nyu", db="scop")
    cursor = db.cursor()
    
    # execute SQL statement
    cursor.execute("SELECT id, pdb_id, chain_id from pdb2scop")
    # pdb_id, chain_id, 
    # (pdb_id, chain_id, sf_start, sf_end, sccs, length, ascession_date)
    all = cursor.fetchall()
    errors=0
    incorrect=0
    for record in all:
        try:
            id, pdb_id, chain_id = record
            
            if len(chain_id.split(":"))>1:
                chain, range = chain_id.split(":")
            else:
                chain, range = (chain_id, "")
            start, end = (None, None)
            if len(range.strip()) > 0:
                try:
                    left, right = range.strip().replace("-"," ").split()
                except:
                    raise Exception("Range no good range %s pdb_id %s chain_id %s"%(range, pdb_id, chain_id))
                start, end = (int(left), int(right))
            
            seq_id = "%s_%s" % (pdb_id, chain)
            cursor.execute("SELECT s.sequence FROM sequence s where s.pdb_id='%s'" % (seq_id))
            record = cursor.fetchone()
            if record:
                sequence = record[0]
            else:
                #Try without trailing Chain letter
                seq_id = "%s_" % pdb_id
                cursor.execute("SELECT s.sequence FROM sequence s where s.pdb_id='%s'" % (seq_id))
                record = cursor.fetchone()
                if record:
                    sequence = record[0]
                else:            
                    raise Exception("Record null for pdb_id %s and chain_id %s and seq_id %s" % (pdb_id, chain_id, seq_id))
            if start and end:
                #print "Trimming from %i to %i from length %i" % (start, end, len(domain_sequence))
                domain_sequence = sequence[abs(start):abs(end)]
            else:
                domain_sequence = sequence
            
            if not len(domain_sequence) > 0:
                if len(sequence) > 0:
                    incorrect +=1
                    raise Exception("Parsing is bad ",start,end,pdb_id, chain_id, seq_id, sequence, domain_sequence)
                else:
                    raise Exception("Domain Sequence is empty ",start,end,pdb_id, chain_id, seq_id, sequence, domain_sequence)
            
            #print ">%s|%s|%s\n%s" % (id, pdb_id, chain_id, domain_sequence)
        except ValueError:
            print "Chain", chain_id
        except Exception, e:
            errors+=1
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
            #raise
    print >>sys.stderr, "Errors ",errors,"Incorrect",incorrect
    