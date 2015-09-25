import os
import subprocess
from Bio.SeqRecord import SeqRecord
from hpf.utilities import QueueParser

CREATE_TABLE = """
 CREATE TABLE %s (
  `id` int(11) NOT NULL auto_increment,
  `cluster_id` int(11) NOT NULL default '0',
  `cluster_entry_num` int(11) NOT NULL default '0',
  `aa_count` int(11) NOT NULL,
  `percent` double unsigned NOT NULL,
  `identifier` varchar(10) default NULL,
  `representative` bool default FALSE,
  `timestamp` timestamp NOT NULL default NOW(),
  `experiment` varchar(20) default NULL,
  `foreign_key` int(11) DEFAULT NULL,
  PRIMARY KEY  (`id`),
  KEY `cluster_id` (`cluster_id`),
  KEY `identifier` USING BTREE (`identifier`),
  KEY `experiment` (`experiment`),
  KEY `foreign_key` (`foreign_key`)
) 
"""

INSERT_ENTRY = """
(cluster_id, cluster_entry_num, aa_count, percent, identifier, representative, experiment, foreign_key)
values
(%s,%s,%s,%s,%s,%s,%s,%s)
"""

def insert_format(entry, experiment=None, foreign_key=None):
    """Return a tuple formatted for the insert statement"""
    return (entry.cluster_id,entry.num,entry.aa,entry.perc,entry.id,entry.repr, experiment, foreign_key)

#ENGINE=MyISAM AUTO_INCREMENT=1104575 DEFAULT CHARSET=latin1

class CDHitOptions(object):
    
    def __init__(self, input, output=None, identity=0.8, length=0.8):
        self._itemp=None
        self._otemp=None
        self.identity=identity
        self.length=length
        if isinstance(input, str):
            assert(os.path.exists(input))
            self.input = input
        elif isinstance(input, list):
            from tempfile import NamedTemporaryFile
            from Bio.SeqIO.FastaIO import FastaWriter
            self._itemp = NamedTemporaryFile()
            self.input = self._itemp.name
            writer = FastaWriter(self._itemp,wrap=0)
            writer.write_records(input)
            self._itemp.flush()
        else:
            raise Exception("Unknown input type",input)
        
        if isinstance(output, str):
            self.output = output
        elif output==None:
            self._otemp = NamedTemporaryFile()
            self.output = self
            
    def __del__(self):
        for temp in [self._itemp, self._otemp]:
            if temp:
                temp.close()

class CDHit(object):
    
    def __init__(self, options):
        self.options = options
        
    def run(self):
        cmd = "cd-hit -i %s -o %s -c %f" % (self.options.input, self.options.output, self.options.identity)
        if self.options.length:
            cmd += " -s %f" % (self.options.length)
        with open(os.devnull) as handle:
            subprocess.check_call(cmd, shell=True, stdout=handle, stderr=handle)
        return CDHitParser(self.options.output+".clstr")
    
class CDHitCluster(list):
    
    def __init__(self, id):
        self.id = id
        
    def representative(self):
        for r in self:
            if r.repr:
                return r
        return None

class CDHitEntry(object):
    
    def __init__(self,id,aa,num,perc,cluster_id,repr=False):
        self.id = id
        self.aa = aa
        self.num = num
        self.perc = perc
        self.cluster_id = cluster_id
        self.repr = repr

class CDHitParser(QueueParser):
    """Parses cd-hit log files"""

    def next(self):
        #>Cluster 16
        #0       234aa, >8546205... at 100%
        #1       331aa, >8502010... *
        if(self._finished()):
            raise StopIteration()
        line = self._get()
        if not line.startswith(">"):
            raise Exception("Cluster",line)
        cluster = CDHitCluster(int(line.split()[-1]))
        while True:
            if self._peek().startswith(">") or self._finished():
                break
            line = self._get().replace("\t"," ")
            parts = line.split()
            num = int(parts[0])
            aa = int(parts[1][:-3])
            id = int(parts[2][1:-3])
            if len(parts) == 5:
                perc = int(parts[4][:-1])
                repr = False
            elif len(parts) == 4:
                perc = 100
                repr = True
            else:
                raise Exception("Entry",line)
            cluster.append(CDHitEntry(id,aa,num,perc,cluster.id,repr))
        if len(cluster) == 1:
            cluster[0].repr=True
        return cluster
    
