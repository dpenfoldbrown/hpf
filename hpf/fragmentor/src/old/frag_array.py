#!/usr/bin/python
import sys, getopt, os, time
from hpf import utilities as util
from hpf.utilities import paths
import fileinput

# Home and scripts
pl_script = 'make_fragments.local.pl'
home = ""
shareware = ""
results=""

task_id = None
code = ""
name = ""
fasta_file = ""

# Scratch location
scr = '/scratch/fragmentor'
scr_job = ""
scr_fasta = ""
cleanup = True
no_homs = None

def usage():
    print "usage: ", sys.argv[0], "code fasta_file"
    print "code = The two letter code"
    print "fasta_file = The fasta file."
    print "results_dir = The results dir minus code and name and all that"
    print "              defaults to RESULTS_DIR/[code]/[code][task_id]"
    print ""
    print "OPTIONS"
    print "-t = the taskId, otherwise taken from PBS_ARRAYID-1"
    print "-n = specifies -nohoms to fragment picker."
    print "-p = The properties file to use for specifying things like the scratch dir, results path."
    print "-c = DO NOT cleanup, defaults to true"
    print ""
    print "NOTES ON FILES"
    print "Inject the task id into path at %s"
    print "  example, -t 0 -o /tmp/lj%s.fasta is replaced with /tmp/lj000.fasta"
    print ""
    print "NOTES"
    print "Assumes nr and nnmake data is in the scratch directory"

def replicate_fasta():
    global scr_fasta
    
    file = paths.getFile(fasta_file)
    scr_fasta = paths.join(scr_job, file)
    if not paths.exists(scr_fasta):
        util.debug("Replicating fasta file")
        util.copy(fasta_file, scr_fasta)
        
    paths.existsOrFail(scr_fasta)

def options(argv):
    global home, shareware,scr, scr_job, fasta_file, code, results, name, cleanup, task_id, no_homs

    if len(argv) < 3:
        usage()
        sys.exit(2)
        
    home = os.getcwd()
    shareware = paths.join(home, "tools")
    
    opts, args = getopt.getopt(argv, "t:p:cn")
    code = args[0]
    fasta_file = args[1]
    results = args[2]
    
    for o,a in opts:
        if o == '-t':
            id ='%(#)03d' % {"#": int(a)}
            task_id = id
        if o == '-p':
            loadProps(a)
        if o == '-c':
            cleanup = False
        if o == '-n':
            no_homs = True
    
    if task_id == None:
        try:
            id ='%(#)03d' % {"#": int(os.environ['PBS_ARRAYID'])-1}
            print "ID",id
            task_id = id
        except:
            util.error("NO TASK ID")

    name = "%s%s" % (code, task_id)
    print code," ",task_id," ",name
    try:
        fasta_file = fasta_file.replace("%s", name)
    except:
        raise
    
    results = paths.join(results, code, name)
    scr_job = paths.join(scr, "jobs", name)
    return

def loadProps(file):
    pass

def runScript():
    # Script requires the following environment variables
    # BLAST_DIR
    # NR_DIR
    # NNMAKE_DIR - the BLOSUM score matrices
    # PSIPRED_DIR
    # NNMAKE_DIR
    # JUFO_DIR
    # SAM_DIR
    # SAM_2ND_DIR - SAM Secondary Structure Prediction
    util.debug("Running fragmentation script")
    
    env = {
           'BLAST_DIR':'%s/blast' % shareware,
           'NR_DIR':'%s/db/nr' % scr,
           'NNMAKEDB_DIR':'%s/db/nnmake_database' % scr,
           'NNMAKE_SHORT_DIR':'%s/nnmake' % shareware,
           'PSIPRED_DIR':'%s/psipred' % shareware,
           'JUFO_DIR':'%s/jufo' % shareware,
           'SAM_DIR':'%s/sam' % shareware,
           'SAM_2ND_DIR':'%s/sam.predict-2nd' % shareware
           }
    
    for key in env :
        # Make sure all paths in 'env' exist in the system.
        paths.existsOrFail(env[key])
        # Put 'env' values into the system environment.
        os.environ[key] = env[key] 
    
    script = paths.join(home, "scripts", pl_script)
    if no_homs:
        script = "%s -nohoms" % script
    fasta = paths.getFile(scr_fasta)
    cmd = "cd %s && %s -verbose -nosam %s" % (scr_job, script, fasta)
    #cmd = "cd %s && %s -verbose %s" % (scr_job, script, fasta)
    util.system(cmd)

def replicate():
    paths.ensure(scr_job)
    replicate_fasta()

def gather():
    util.debug("Gathering result")
    job_results = results
    paths.ensure(job_results)
    
    
    #aazj00103_05.075_v1_3.bz2  zj001.fasta    zj001.psipred_ss2
    #aazj00109_05.075_v1_3.bz2  zj001.psipred
    aa_pattern = "aa"+name+"(03|09)\_05\.075\_v1\_3$"
    patterns = [aa_pattern]
    for suf in [".fasta", ".psipred", ".psipred_ss2"]:
        patterns.append(name+suf)
    
    try:
        for pattern in patterns:
            print "Pattern ",pattern
            match = paths.find(pattern, scr_job)
            if not match:
                raise "Missing file"
            for path in match:
                if pattern == aa_pattern:
                    print "Filtering Columns %s" % path
                    for line in fileinput.input(path, inplace=1):
                        print " "+line[:47].strip()
                    fileinput.close()
                print "Path ", path
                file = paths.getFile(path)
                dest = paths.join(job_results, file)
                util.copy(path, job_results)
                if pattern == aa_pattern:
                    util.system("bzip2 %s"  % dest)
    except:
        paths.removerf(job_results)
        raise
        
def cleanup():
    paths.removerf(scr_job)

def main(argv):
    try:
        options(argv)
    except:
        usage()
        raise
        sys.exit(2)
    
    try:
        replicate()
        runScript()
        gather()
    finally:
        if cleanup:
            cleanup()
    
    util.debug("FINISHED!")
    
if __name__ == "__main__":
    print "YAY"
    main(sys.argv[1:])
