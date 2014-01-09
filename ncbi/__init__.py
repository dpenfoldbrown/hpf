from Bio import Entrez

MAX_RETURN=100000

def gis(taxonomy_id):
    count = count_gi(taxonomy_id)
    start = 0
    while start<count:
        handle = Entrez.esearch(db="protein",term="txid%i[Organism:exp]" % taxonomy_id,retstart=start,retmax=MAX_RETURN)
        try:
            record = Entrez.read(handle)
            for id in record['IdList']:
                yield int(id)
            start+=MAX_RETURN
        finally:
            handle.close()

def count_gi(taxonomy_id):
    """Return the number of GI's for an organism id"""
    handle = Entrez.esearch(db="protein",term="txid%i[Organism:exp]" % taxonomy_id,retstart=0,retmax=2)
    try: 
        record = Entrez.read(handle)
        count = int(record['Count'])
        return count
    finally:
        handle.close()

class BlastPGPTask(OutputTask):
    """
    Blast a given input to the HPF database
    """
    
    def __init__(self, input="input.aln", db="nr", output="blast.xml", **kwargs):
        OutputTask.__init__(self, *output)
        self.input = intput
        self.db = db
        self.kwargs = kwargs

    def _do(self):
        import subprocess
        blast_exe = subprocess.Popen("which blastpgp", shell=True, stdout=subprocess.PIPE).communicate()[0].strip()
        runtime().debug("Blasting with alignment %s using %s" %(self.input,blast_exe))
        r,e = NCBIStandalone.blastpgp(blast_exe, 
                                      self.db,
                                      self.alignment,
                                      **self.kwargs)
        consume(r)
        return tuple(self.output)