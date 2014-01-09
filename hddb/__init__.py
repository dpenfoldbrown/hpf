import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from hpf.runtime import runtime
from Bio import SeqIO

SCRIPTS_FOLDER = os.path.join(os.path.dirname(__file__),"scripts")

def tunnel(sleep=None):
    """
    Ensure a tunnel to the uwashinton mysql server is running
    """
    from scripts import TUNNEL
    if sleep!=None:
        import time, random
        time.sleep(random.random()*sleep)
    runtime().debug(TUNNEL)
    import subprocess
    with open(os.devnull) as handle:
        subprocess.call(TUNNEL,shell=True,stdout=handle,stderr=handle)

def kill_tunnel():
	import subprocess
	from scripts import KILLTUNNEL
	with open(os.devnull) as handle:
		subprocess.call(KILLTUNNEL,shell=True,stdout=handle,stderr=handle)

def connect(db="bddb", user="pwinters", passwd="bonneaulab", **kwargs):
    """
    Ensures a tunnel is connected, then returns a MySQLdb connection object to the hddb.
    """
    tunnel()
    import MySQLdb
    return MySQLdb.connect(host="127.0.0.1",port=13307,user=user,passwd=passwd,db=db,**kwargs)

def ginzu_svg(sequence_key, width=800):
    """
    Return the vector graphic for GINZU domains.
    """
    import subprocess
    cmd = "perl %s %i %i" % (os.path.join(SCRIPTS_FOLDER,"svg.pl"), sequence_key, width)
    runtime().debug(cmd)
    print cmd
    svg = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]+"\n"
    return svg

def sequences(cursor, sequence_key=None):
    """
    Return the given sequences as SeqRecord objects.
    """
    query = "SELECT id,sequence from sequence"
    if sequence_key != None:
        if not hasattr(sequence_key, "__iter__"):
            sequence_key = [sequence_key]
        query += " where id in (%s)" % (",".join([str(key) for key in sequence_key]))
    runtime().debug(query)
    cursor.execute(query)
    runtime().debug(query)
    for id, sequence in cursor.fetchall():
        yield SeqRecord(Seq(sequence), str(id))
        
def proteins(cursor, experiment=None, filter_experiments=True, sequence_key=None):
    """
    Return the selected proteins as SeqRecord objects
    """
    query = """SELECT s.id,s.sequence, e.id, e.short_name, e.taxonomy_id
        from hpf.experiment e 
        join bddb.protein p on e.id=p.experiment_key
        join ddbCommon.sequence s on p.sequence_key=s.id 
        """
    assert experiment!= None or sequence_key != None
    if experiment != None or filter_experiments==True or sequence_key != None:
        query += " where "
    if experiment:
        if not hasattr(experiment, "__iter__"):
            experiment = [experiment]
        query += " e.id in (%s)" % (",".join([str(key) for key in experiment]))
    if filter_experiments:
        t = " e.taxonomy_id!=0"
        query += " and "+t if experiment else t
    if sequence_key:
        t = " s.id in (%s)" % (",".join([str(key) for key in sequence_key]))
        query += " and "+t if experiment or filter_experiments else t
    runtime().debug(query)
    cursor.execute(query)
    runtime().debug("Fetching")
    for id, sequence, e_id, e_name, taxonomy_id in cursor.fetchall():
        record = SeqRecord(Seq(sequence), str(id), description=e_name)
        record.annotations = {"taxonomy_id":taxonomy_id,
                              "experiment_key":e_id,
                              "organism":e_name}
        yield record
    
