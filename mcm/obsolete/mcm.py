#!/usr/bin/env python

import sys,os
import getopt
import time
import subprocess
import shutil
import re
from hpf.runtime import runtime, Flag, IntegerOption, FileOption
from hpf.utilities import consume
from hpf.processing import processor, SYNCHRONOUS
from hpf.hddb.mcm import FOLDER,SCRIPTS_FOLDER,TOOLS_FOLDER

SCRATCH = "scratch"
RESULTS = "results"
DATABASE = "database"
TRIM_RG = "trim_rg"
TRIM_SCORE = "trim_score"
CLEANUP = "cleanup"

PL_SCRIPT = os.path.join(SCRIPTS_FOLDER,"mcm.pl")

def _mcm(args):
    print args
    return mcm(*args)

def mcm(decoy_file, psipred_file):
    m = MCM(decoy_file,
        psipred_file,
        pred_code(decoy_file),
        runtime().opt(SCRATCH),
        results = runtime().opt(RESULTS),
        database = runtime().opt(DATABASE),
        trim_score = runtime().opt(TRIM_SCORE),
        trim_rg = runtime().opt(TRIM_RG),
        cleanup = not runtime().opt(CLEANUP))

    m.run()
    
class MCM(object):
    
    def __init__(self,
                 decoy_file,
                 psipred_file,
                 job_name,
                 scratch,
                 results=None,
                 database=None,
                 trim_score=15000,
                 trim_rg=10000,
                 cleanup=False):
        self.decoy_file = decoy_file
        self.psipred_file = psipred_file
        self.job_name = job_name
        self.scratch = scratch
        self.results = results if results else scratch
        self.database = database if database else scratch
        self.job = os.path.join(scratch,job_name)
        self.report_file = None
        self.trim_score = trim_score
        self.trim_rg = trim_rg
        self.cleanup = cleanup
    
    def run(self):
        try:
            prep_time = time.localtime()
            self.report("Preparation:%s" % time.asctime(prep_time))
            self._replicate()
            if self.trim_score or self.trim_rg:
                self._trim()
            self.start_time = time.localtime()
            self.report("Begin:%s" % time.asctime(self.start_time))
            self._runScript()
            self.end_time = time.localtime()
            self.report("End:%s" % time.asctime(self.end_time))
            self.report("Total Seconds:%d" % (time.mktime(self.end_time)-time.mktime(self.start_time)))
            self._gather()
        finally:
            if self.cleanup:
                self._cleanup()
        
        
    def _trim(self):
        global scr_outfile
        
        scriptsDir = SCRIPTS_FOLDER
        file = self.decoy_file
        trimmedFile = self.decoy_file+".trim" 
        renameFile = trimmedFile[:-len(".out.trim")]
        renameFile += ".%d.%d.out" % (self.trim_score, self.trim_rg)
        
        if not os.path.exists(renameFile) :
            score = "-s %d " % self.trim_score if self.trim_score != 0 else ""
            rg = "-r %d" % self.trim_rg if self.trim_rg != 0 else ""
            cmd = "%s/trim_score_rg.pl -f %s %s %s" % (scriptsDir, file, self.trim_score, self.trim_rg)
            subprocess.check_call(cmd,shell=True)
            os.rename(trimmedFile, renameFile)
            
        self.decoy_file = renameFile
    
    def _runScript(self):
        # Script variables
        cmd = "%s -outfile %s -ssfile %s -epath %s -dpath %s" % (PL_SCRIPT, 
                                                                 self.decoy_file, 
                                                                 self.psipred_file, 
                                                                 TOOLS_FOLDER, 
                                                                 self.database)
        subprocess.check_call(cmd,cwd=self.job,shell=True)
    
    def _extract(self, file):
        """Extract an input file and filter with dos2unix"""
        from hpf.utilities.paths import isZipped, zipSuffix
        from hpf.utilities import extract
        #
        destination = os.path.join(self.job,os.path.basename(file))
        #extract = paths.isZipped(file)
        #rename = not file.endswith(".out")
    
        if not os.path.exists(destination) :
            shutil.copy(file, destination)
            extracted = ".".join(destination.split(".")[:-1]) if isZipped(destination) else destination
            if isZipped(destination) and not os.path.exists(extracted):
                destination, suffix = extract(destination)
        cmd = "dos2unix %s" % destination
        subprocess.check_call(cmd, shell=True, cwd=self.job)
        assert os.path.exists(destination)
        return destination
        
    def _replicate(self):
        """Copy input files to the job directory ins scratch"""
        global scr_job
        from hpf.utilities.paths import ensure
        ensure(self.job)
        self.decoy_file = self._extract(self.decoy_file)
        if not self.decoy_file.endswith(".out"):
            destination = self.decoy_file+".out"
            shutil.move(self.decoy_file, destination)
            self.decoy_file = destination
        self.psipred_file = self._extract(self.psipred_file)
        
    def _gather(self):
        """Gather result files and report text file."""
        from hpf.utilities.paths import ensure, find
        scr_results_dir = self.decoy_file+".p"
        results = os.path.join(self.results,self.job_name)
        ensure(results)
        # Using the full filename will ensure the directory is created.
        shutil.copy(os.path.join(scr_results_dir, "log.xml"), results)
        
        decoys = list(find("decoy\_[0-9]+\.pdb",dir=scr_results_dir))
        assert len(decoys)>0
        for decoy in decoys:
            shutil.copy(decoy, results)
        shutil.copy(os.path.join(self.job, "report.txt"), results)
        
    def report(self, string):
        """Append a string to the report file"""
        if not self.report_file:
            from hpf.utilities.paths import ensure
            ensure(self.job)
            self.report_file = open(os.path.join(self.job,"report.txt"), "a")
        print >>self.report_file, string
        
    def _cleanup(self):
        shutil.rmtree(self.job)
        #paths.removerf(self.job)

def pred_code(decoy):
    r = re.compile("[a-z]{2}[0-9]{3}\.((result)|(out))\.gz")
    match = r.search(decoy)
    return decoy[match.start():match.end()][:5]

def tasks(decoy_dir, psipred_dir):
    from hpf.utilities.paths import find
    runtime().debug("Searching for decoys",decoy_dir)
    decoys = list(find("[a-z]{2}[0-9]{3}\.((result)|(out))\.gz",dir=decoy_dir))
    runtime().debug("Found",len(decoys),"decoys")
    for decoy in decoys:
        prediction_code = pred_code(decoy)
        psipred = os.path.join(psipred_dir,prediction_code,prediction_code+".psipred")
        if os.path.exists(psipred):
            yield (decoy,psipred)
        

def main(*args):
    pool = processor()
    runtime().debug("Using processor",pool)
    pool.make_tasks(tasks,*args)
    consume(pool.run(_mcm))
    
def _do(argv):
    r = runtime()
    r.description("""
    mcm.py [-options] decoy_dir, psipred_dir
    template script.
    """)
    r.add_option(Flag(CLEANUP, "c", description="Leave scratch directory on error.", default=False))
    r.add_option(IntegerOption(TRIM_SCORE,"s",description="Sort decoys by score and keep #", default=15000))
    r.add_option(IntegerOption(TRIM_RG,"g",description="Sort decoys by radius of gyration and keep # after score filter", default=10000))
    r.add_option(FileOption(SCRATCH,"x",description="Scratch/job directory to use.",default="/scratch/pcw216/tmp/mcm/jobs"))
    r.add_option(FileOption(RESULTS,"r",description="Results directory to move files to.",default="/scratch/pcw216/tmp/mcm/results"))
    r.add_option(FileOption(DATABASE,"z",description="Database directory to use.",default="/scratch/pcw216/tmp/mcm/db"))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
