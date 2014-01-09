#!/usr/bin/env python
'''
Created on Apr 6, 2010

Takes a list of prediction codes to set up input files and create a task pickle 
for running MCM.
Originally used to re-do the since solved.

@author: Patrick
'''
from hpf.hddb import tunnel
tunnel(sleep=1)
import sys, os
import getopt
import shutil
from hpf.runtime import runtime, Flag, FileOption, debug, DEBUG
from hpf.utilities import consume
from hpf.processing import processor, SYNCHRONOUS
from hpf.mcm.mcm import SCRATCH
from hpf.hddb.db import *

FRAGMENTS = "fragments"
SCRATCH = "SCRATCH"
HPF1_DECOYS = "hpf1_decoys"
HPF2_DECOYS = "hpf2_decoys"

def prep(prediction_code):
    runtime().debug("Running on",prediction_code)
    with PrepPredictionCode(prediction_code,
                            runtime().opt(SCRATCH),
                            runtime().opt(HPF1_DECOYS),
                            runtime().opt(HPF2_DECOYS)) as p:
        return p()

class PrepPredictionCode(object):
    
    def __init__(self, 
                 prediction_code,
                 scratch,
                 hpf1=None,
                 hpf2=None):
        """
        @param hpf#: The path to decoys for HPF1 and HPF2 respectively.
            A find and replace will be put on the path to define variables:
                %prediction_code%
                %prediction_letter%
            ie hpf2="/archive/rb133/hpf2/%prediction_letter%_hpf2_results/%prediction_code%.result.gz"
        """
        self.prediction_code = prediction_code
        self.scratch = os.path.abspath(scratch)
        self.hpf1 = hpf1
        self.hpf2 = hpf2
    
    def __enter__(self):
        from hpf.utilities.paths import ensure
        ensure(self.scratch)
        self.session = Session()
        runtime().debug("Loading outfile")
        self.filesystem_outfile = self.session.query(FilesystemOutfile).filter(FilesystemOutfile.prediction_code==self.prediction_code).first()
        runtime().debug("Loading sequence")
        self.sequence = self.session.query(Sequence).get(self.filesystem_outfile.sequence_key)
        debug(self.prediction_code,self.filesystem_outfile,self.sequence)
        return self
        
    def __exit__(self, type, value, traceback):
        self.session.close()
    
    def __call__(self, *args, **kwargs):
        return self.run()
    
    def run(self):
        from hpf.utilities.paths import ensure
        
        if self.filesystem_outfile.executable_key==179:
            runtime().debug("HPF1 format",self.filesystem_outfile)
            _format = self.hpf1
        else:
            runtime().debug("HPF2 format",self.filesystem_outfile)
            _format = self.hpf2

        debug("Using decoy format",_format)
        self.decoy = format(_format, self.prediction_code)
        if not os.path.exists(self.decoy):
            debug("NO DECOY FILE",self.prediction_code,self.decoy)
            #return None
        else:
            self.destination = os.path.join(self.scratch,os.path.basename(self.decoy))
            if not os.path.exists(self.destination):
                ensure(os.path.join(self.scratch,self.prediction_code))
                self.decoy = self._decoy()
            else:
                debug("exists",self.destination)
                self.decoy = self.destination
        
        self.psipred = os.path.join(self.scratch,self.prediction_code,self.prediction_code+".psipred")
        if not os.path.exists(self.psipred):
            ensure(os.path.join(self.scratch,self.prediction_code))
            self.psipred = self._ss()
        else:
            debug("exists",self.psipred)

        if all(map(os.path.exists,[self.decoy,self.psipred])):
            runtime().debug("Exported",(self.decoy,self.psipred))
            return (self.decoy,self.psipred)
        else:
            runtime().debug("Failed to export")
            return None
        
    def _decoy(self):
        """
        Find the decoy file and copy it to the scratch directory.
        """

        assert os.path.exists(self.decoy), self.decoy
        debug("Found",self.decoy)
        destination = self.destination
        shutil.copy(self.decoy, destination)
        #assert os.path.exists(destination)
        debug("Copied to",destination)
        return destination
    
    def _ss(self):
        self.psipred = self._db() if self.sequence.psipred != None else self._psipred()
        assert os.path.exists(self.psipred)
        return self.psipred
    
    def _db(self):
        """
        Get the prediction from the database and write it to a file.
        """
        debug("Getting psipred from database")
        from hpf.pdb.psipred import PsipredWriter
        psipred = self.psipred
        with open(psipred,"w") as handle:
            PsipredWriter().write(handle, 
                                  self.sequence.psipred.psipred, 
                                  self.sequence.record)
        #assert os.path.exists(psipred)
        return psipred
    
    def _psipred(self):
        """
        Run psipred on the sequence.
        """
        debug("Running psipred prediction")
        from hpf.seq import TemporaryRecordFile
        from hpf.pdb.psipred import Psipred, PsipredOptions
        from Bio import SeqIO
        fasta = os.path.join(self.scratch,self.prediction_code,self.prediction_code+".fasta")
        with open(fasta,"w") as handle:
            SeqIO.write([self.sequence.record], handle, "fasta")
        psipred = self.psipred
        from Bio.Blast.NCBIStandalone import blastpgp
        chk = fasta+".chk"
        import subprocess
        cmd = subprocess.Popen(["which", "blastpgp"], stdout=subprocess.PIPE).communicate()[0].strip()
        debug("Using",cmd)
        result,error = blastpgp(cmd, 
                                "nr", 
                                fasta, 
                                npasses=3,
                                checkpoint_outfile=chk,
                                expectation=1e-4,
                                model_threshold=1e-4,
                                align_outfile="/dev/null")    
        debug(result.readlines())
        debug(error.readlines())

        options = PsipredOptions(fasta,
                                 profile=chk,
                                 output=psipred+".1",
                                 output2=psipred+".2",
                                 horiz=psipred,
                                 cwd = os.path.join(self.scratch,self.prediction_code))
        prediction = Psipred(options).run()
        db_pred = PsipredFactory().create(prediction,sequence_key=self.sequence.id)
        self.session.add(db_pred)
        self.session.commit()
        #assert os.path.exists(psipred)
        return psipred

def tasks(*args):
    def codes(input):
        for line in input:
            yield line.split()[0].strip()
    
    if len(args)==0:
        return codes(sys.stdin)
    elif len(args)==1 and os.path.exists(args[0]):
        with open(args[0]) as input:
            return codes(input)
    else:
        return list(args)

def format(decoy_file, prediction_code):
    file = decoy_file
    file = file.replace("%prediction_code%", prediction_code)
    file = file.replace("%prediction_letter%", prediction_code[0:2])
    return file
    
def main(*args):
    runtime().set_debug(1)
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS), raise_errors=False)
    runtime().debug("Using processor",pool)
    pool.make_tasks(tasks)
    mcm_tasks = [t for t in pool.run(prep) if t !=None]
    import cPickle
    with open("exported.pickle","w") as handle:
        cPickle.dump(mcm_tasks,handle)
    
    
def _do(argv):
    r = runtime()
    r.description("""
    template.py [-options] args
    template script.
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing"))
    r.add_option(FileOption(SCRATCH, "x", description="Scratch output directory.", default="exported"))
    r.add_option(FileOption(FRAGMENTS, "f", description="Location of HPF2 fragments", default="/archive/rb133/hpf2/fragments"))
    r.add_option(FileOption(HPF2_DECOYS, "2", description="Location of HPF2 decoys", default="/archive/rb133/hpf2/results/%prediction_letter%_hpf2_results/%prediction_code%.result.gz"))
    r.add_option(FileOption(HPF1_DECOYS, "1", description="Location of HPF1 decoys", default="/archive/rb133/hpf1/results/%prediction_letter%/%prediction_code%.out.gz"))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
