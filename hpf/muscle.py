import os
import subprocess
from Bio.SeqRecord import SeqRecord
from hpf.utilities import QueueParser,OutputTask
from hpf.runtime import runtime

class MuscleOptions(object):
    
#    -in <inputfile>    Input file in FASTA format (default stdin)
#    -out <outputfile>  Output alignment in FASTA format (default stdout)
#    -diags             Find diagonals (faster for similar sequences)
#    -maxiters <n>      Maximum number of iterations (integer, default 16)
#    -maxhours <h>      Maximum time to iterate in hours (default no limit)
#    -maxmb <m>         Maximum memory to allocate in Mb (default 80% of RAM)
#    -html              Write output in HTML format (default FASTA)
#    -msf               Write output in GCG MSF format (default FASTA)
#    -clw               Write output in CLUSTALW format (default FASTA)
#    -clwstrict         As -clw, with 'CLUSTAL W (1.81)' header
#    -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)
#    -quiet             Do not write progress messages to stderr
#    -stable            Output sequences in input order (default is -group)
#    -group             Group sequences by similarity (this is the default)
#    -version           Display version information and exit

    
    def __init__(self, input, output=None, diags=None, maxhours=None, maxmb=None, clwstrict=False, quiet=True, clw=False):
        self._itemp=None
        self._otemp=None
        self.maxhours=maxhours
        self.clwstrict=clwstrict
        self.clw=clw
        self.quiet=quiet
        self.maxmb = maxmb
        if isinstance(input, str):
            assert(os.path.exists(input)),"Can't find file %s" % input
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
            self.output = self._otemp.name
            
    def __del__(self):
        for temp in [self._itemp, self._otemp]:
            if temp:
                temp.close()
    

class Muscle(object):
    
    def __init__(self, options):
        self.options = options
        
    def run(self):
        cmd = "muscle -in %s -out %s" % (self.options.input, self.options.output)
        if self.options.clwstrict:
            cmd += " -clwstrict"
        if self.options.clw:
            cmd += " -clw"
        if self.options.quiet:
            cmd += " -quiet"
        if self.options.maxhours:
            cmd += " -maxhours %i" % self.options.maxhours
        if self.options.maxmb:
            cmd += " -maxmb %i" % self.options.maxmb
        subprocess.check_call(cmd, shell=True)
        return self.options.output
    
class MuscleTask(OutputTask):
    """
    Perform Muscle multiple sequence alignment.  Limits to running for 1 hour.
    """

    def __init__(self,fasta,output,phylip=False,**kwargs):
        OutputTask.__init__(self,output)
        self.fasta = fasta
        self.phylip = phylip
        self.kwargs = kwargs

    def _do(self):
        muscle = Muscle(MuscleOptions(input=self.fasta,output=self.output,maxhours=1,**self.kwargs))
        runtime().debug("Performing Muscle w/", self.fasta)
        out = muscle.run()
        if self.phylip:
            from hpf.phylip import clustal_to_phylip
            return clustal_to_phylip(self.output)
        else:
            return self.output
        
