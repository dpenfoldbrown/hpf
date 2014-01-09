#!/usr/bin/env python

# WARNING NOT REALLY WORKING

# This script will process, parse, and store in the database all of the ginzu cuts for domain information
#CUTS      q_beg q_end q_len m_beg m_end m_len p_beg p_end    p_id        conf     source m_seq
#  1   116   116     1   118   118     0     0      na   52.357997        msa GSHMTALDWRSALTADEQRSVRALVTATTAVDGVAPVGEQVLRELGQQRTEHLLVAGSRPGGPIIGYLNLSPPRGAGGAMAELVVHPQSRRRGIGTAMARAALAKTAGRNQFWAHGTL
#117   318   202   115   318   204     2   174  1ufhA_    3.154902   pdbblast HGTLDPARATASALGLVGVRELIQMRRPLRDIPEPTIPDGVVIRTYAGTSDDAELLRVNNAAFAGHPEQGGWTAVQLAERRGEAWFDPDGLILAFGDSPRERPGRLLGFHWTKVHPDHPGLGEVYVLGVDPAAQRRGLGQMLTSIGIVSLARRLGGRKTLDPAVEPAVLLYVESDNVAAVRTYQSLGFTTYSVDTAYALAGTDN


import MySQLdb
import getopt
import sys
from collections import defaultdict
from hpf.processing import MultiProcessor

db = None

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "upload"
opts = {DEBUG: False, SILENT: False, UPLOAD: True}

def cuts(args):
    """Use the ginzu key and, given domains, parse the cuts file"""
    protein_key, domains = args
    
    #print "RUN CUTS", args
    
    global db
    if db==None:
        db = _connection()
    values = []
    cursor = db.cursor()
    try:
        cursor.execute("select id, uncompress(compress_cuts) from ddbCommon.ginzuRun where sequence_key=%i order by version desc" % protein_key)
        ginzu_key, compressed = cursor.fetchone()
        domain_cuts = compressed.strip().split("\n")[1:]
        # Parse out values
        for line in domain_cuts:
            #print "LINE",type(line),line
            q_beg, q_end, q_len, m_beg, m_end, m_len, p_beg, p_end, p_id, conf, source, m_seq = line.split()
            if p_id == "na":
                p_id = "NULL"
            domain_sequence_key = None
            # Find the appropriate sequence key
            for sequence, dkey in domains:
    #            print sequence
    #            print m_seq
    #            print m_seq.find(sequence) != -1 or sequence.find(m_seq) != -1
    #            print ""
                start = int(q_beg)-int(m_beg)
                #end = int(q_len)
                end = len(m_seq)-abs(int(q_end)-int(m_end))
                if m_seq[start:end] == sequence:
                #if m_seq.find(sequence) != -1 or sequence.find(m_seq) != -1:
                    if domain_sequence_key != None and dkey != domain_sequence_key:
                        #print domain_cuts
                        raise Exception("Domain already assigned", m_seq,domain_sequence_key,sequence,dkey, ginzu_key, domains)
                    domain_sequence_key = dkey
                
            
            if domain_sequence_key == None:
                dm = m_seq[start:end]
                match = False
                for sequence,dkey in domains:
                    if (sequence.find(m_seq[start:end]) != -1 or m_seq[start:end].find(sequence) != -1):
                        match = True
                        print ""
                        print q_beg, q_end, q_len, m_beg, m_end, m_len
                        print start,end,len(m_seq)
                        print ginzu_key
                        print ""
                        print sequence, len(sequence)
                        print m_seq, len(m_seq)
                        print dm, len(dm)
                        print sequence.find(dm) != -1 or dm.find(sequence) != -1
                        print ""
                if match:
                    raise Exception("Domain sequence key is None", ginzu_key, m_seq, dkey, sequence)
                else:
                    continue
                        
            insert_value = (q_beg, q_end, q_len, m_beg, m_end, m_len, p_beg, p_end, p_id, conf, source, domain_sequence_key, ginzu_key)
            values.append(insert_value)
    finally:
        cursor.close()
    return values

def upload(values):
    
    if opts[UPLOAD]:
        #print "Upload",values
        global db
        if db==None:
            db = _connection()
        cursor = db.cursor()
        try:
            cursor.executemany(
                """INSERT INTO hpf.ginzu_cuts (q_beg, q_end, q_len, m_beg, m_end, m_len, p_beg, p_end, p_id, conf, source, domain_sequence_key, ginzu_key)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)""",
                values)
        finally:
            cursor.close()
    return values

    

def _connection():
    return MySQLdb.connect(host="127.0.0.1", user="pwinters", passwd="bonneaulab", db="hpf", port=3307)    

def main():
    print opts
    db = _connection()
    cursor = db.cursor()

    query = """select p.sequence_key, d.domain_sequence_key, s.sequence from hpf.experiment e join bddb.protein p join ddbCommon.domain d join ddbCommon.ddb_sequence s 
        on e.id=p.experiment_key and p.sequence_key=d.parent_sequence_key and d.domain_sequence_key=s.id where d.domain_type in ('fold_recognition', 'psiblast') limit 50"""
    if not opts[SILENT]:
        print query
    cursor.execute(query)
    proteins = defaultdict(lambda: [])
    done = 0
    for protein_sequence_key, domain_sequence_key, sequence in cursor.fetchall():
        proteins[protein_sequence_key].append((sequence, domain_sequence_key))
        done +=1
        if done % 1000 == 0:
            print "Sorting",done
    # Just make sure we close all this down and disconnect everything before multiprocessing uses shared memory
    cursor.close(); db.close(); db=None
    # Separate each ginzu cuts file into its own task
    tasks = [(key, proteins[key]) for key in proteins]
    print "Tasks set up"
    
    # Set up pool and start running
    outstream = (opts[SILENT] and [sys.stdout] or [None])[0]
    if opts[DEBUG]:
        from hpf.processing import MapProcessor
        pool = MapProcessor()
    else:
        pool = MultiProcessor(completed_out=outstream)
    result = (opts[UPLOAD] and [upload] or [None])[0]
    
    all = []
    exceptions = []
    for values in pool.run(cuts, tasks, result):
        if not isinstance(values, Exception):
            all += values
        else:
            exceptions.append(values)
    print "Completed",len(all)
    print "Exceptions",pool.exceptions
    if opts[DEBUG]:
        for e in exceptions:
            print e
    

def _do(argv):
    try:
        opts, args = _options(argv)
    except Exception, e:
        print e
        _usage()
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "dxq")
    global opts
    for o,a in _opts:
        if o == '-d':
            opts[DEBUG] = True
        if o == '-x':
            opts[UPLOAD] = False
        if o == '-q':
            opts[SILENT] = True
    return opts, args
    
def _usage():
    print "$ ginzu_cuts.py"
    print "Process Ginzu logfiles and upload into db."
    print ""
    print "Options"
    print "-d Debug mode, uses a synchronous processor."
    print "-x DO NOT upload results."
    print "-q Silent mode."
    sys.exit(2)

if __name__=="__main__":
    _do(sys.argv[1:])