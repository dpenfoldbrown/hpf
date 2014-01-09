#!/usr/bin/python

## Fragment picking functionality - a wrapper/driver for a fragment picking perl script (itself a
## wrapper for alignment tools and nnmake fragment picking functionality). 
##
## Auth: Patrick Winters (Original). Duncan Penfold-Brown (Revision), 12/17/2010.
##


# TODO: Remove all unnecessary calls to utilities files (eg: paths.join = os.path.join, exactly).
# TODO: Clean up temporary directory creation for storing inputs and results, etc, in scratch dirs. Pointlessly complex now.
# TODO: Revise globals (most unneccessary).
# TODO: What is the '%s' formatting for fasta_file input. Why necessary? Are the fasta input files
#       for each script run named with their code and unique ID? 
# TODO: Use a config file - eliminate the -p tag for properties file. Put frag_script, scr_dir, etc. in config file.
#       (Also, when make_fragments.local is incorporated into here, add tool locations, etc, to config)
# TODO: Change utilities debug use to logger debug.


## Imports

import os
import sys
import time
import getopt
import fileinput
from optparse import OptionParser

# From utilities, using: debug, error, system, copy. Unneccesary? 
from hpf import utilities
# From paths, using: join, exists, getFile, find, ensure, removerf, existsOrFail. Most are direct wrappers of os.path functions.
from hpf.utilities import paths


## Globals

usage_str = """usage: %prog [-n] [-c] [-t TASK_ID] [-p PROPERTIES_FILE]  CODE  FASTA_FILE  RESULTS_DIR
  
  CODE          The two-letter prediction code (used in naming results files).
  FASTA_FILE    Input file containing sequences in FASTA format. Filename must contain a '%s', to be replaced by [CODE][TASK_ID].
  RESULTS_DIR   The directory for storing results (without code and name, etc). Results are moved to RESULTS_DIR/[CODE]/[CODE][TASK_ID].

  NOTE: The task_id (or ARRAY_ID-1, if PBS or SGE is present) is inserted into the input file path at %s.
        EG: in '-t 001 mm /tmp/%s.fasta /tmp/results', input is /tmp/mm001.fasta
  
  NOTE: Assumes that nr and nnmake data is in the scratch directory.
"""

# Script with fragmenting functionality
frag_script = 'make_fragments.local.pl'

# Scratch locations for storing temporary script files
scr_dir = '/scratch/patrick/tmp/fragmentor'
scr_job_dir = ""
scr_fasta = ""

# Required input and output locations (cmd-line positional arguments)
code = ""
task_id
fasta_file = ""
results_dir=""

# Directories for script use
home = ""
tools_dir = ""


## Option Setting and Parsing.

# Sets the command-line options for a passed-in OptionParser object (via optparse, Python 2.3+)
def set_options(parser):

    # Option attributes: action, type, dest, metavar, default, help. See optparse documentation.
    # Defaults: action=store, type=string, dest=[name of the option] help=none
    # Note: the 'dest' value is how the argument is accessed after parsing. EG: (opts, args) = parser.parse_args(); opts.dest_name
    # EG:
    #parser.add_option("-c", "--cloudConfig", dest="cloud_config",
    #                  action="store"
    #                  type="string"
    #                  metavar="FILE",
    #                  default="Whatever",
    #                  help="Designates some things that are important")
    
    parser.add_option("-n", action="store_true", dest="nohoms", default=False,
                      help="Flag. Specifies -nohoms to the fragment picker.")
    parser.add_option("-c", action="store_false", dest="cleanup", default=True,
                      help="Flag. Do NOT clean up. Defaults to TRUE (by default, cleanup will occur).")
    parser.add_option("-t", dest="task_id", type="int", metavar="TASK_ID", default=None,
                      help="The task_id. If not specified, taken from PBS_ARRAYID when present.")
    parser.add_option("-p", dest="properties_file", metavar="FILE",
                      help="CURRENTLY NOT SUPPORTED. The properties file to use for specifying scratch dir, results path, etc.")

# Parses positional arguments necessary for functionality (CODE, FASTA_INPUT, RESULTS_DIR)
# Takes a list of args and a parser Option object.
def parse_positional_args(args, parser):
    global code, fasta_file, results_dir

    if len(args) < 3:
        print "Error: Incorrect usage."
        parser.print_help()
        sys.exit(1)
    
    # Set input and output globals from positional arguments
    # TODO: Do some verifcation on positional inputs - check for file existence, etc.
    code = args[0]
    fasta_file = args[1]
    results_dir = args[2]
    
    return


## Fragment picking functionality

# Sets the environment and runs the perl fragment making wrapper (make_fragments.local.pl)
def runScript(no_homs=False):
    # Script requires the following environment variables
    # BLAST_DIR
    # NR_DIR
    # NNMAKE_DIR - the BLOSUM score matrices
    # PSIPRED_DIR
    # NNMAKE_DIR
    # JUFO_DIR
    # SAM_DIR
    # SAM_2ND_DIR - SAM Secondary Structure Prediction
    utilities.debug("Running fragmentation script")
    
    env = {
           'BLAST_DIR':'%s/blast' % tools_dir,
           'NR_DIR':'%s/db/nr' % scr_dir,
           'NNMAKEDB_DIR':'%s/db/nnmake_database' % scr_dir,
           'NNMAKE_SHORT_DIR':'%s/nnmake' % tools_dir,
           'PSIPRED_DIR':'%s/psipred' % tools_dir,
           'JUFO_DIR':'%s/jufo' % tools_dir,
           'SAM_DIR':'%s/sam' % tools_dir,
           'SAM_2ND_DIR':'%s/sam.predict-2nd' % tools_dir
           }
    
    for key in env :
        # Check if all paths given in env exist.
        paths.existsOrFail(env[key])
        
        # os.environ is the system environment. New env. vars. can be set by adding to it.
        os.environ[key] = env[key] 
    
    script = paths.join(home, "scripts", frag_script)
    
    if no_homs:
        script = "%s -nohoms" % script
    fasta = paths.getFile(fasta_file)
    cmd = "cd %s && %s -verbose -nosam %s" % (scr_job_dir, script, fasta)
    #cmd = "cd %s && %s -verbose %s" % (scr_job_dir, script, fasta)
    utilities.system(cmd)

# Check for results files. If any are missing, raise an exception. 
def gather():
    global code, task_id

    utilities.debug("Gathering result")
    paths.ensure(results_dir)
    
    # Contents of results directory (3' library, 9' library, fasta, psipred, and psipred_ss2 files). EG:
    #aazj00103_05.075_v1_3.bz2  zj001.fasta    zj001.psipred_ss2
    #aazj00109_05.075_v1_3.bz2  zj001.psipred
    
    # Set "patterns" for results files.
    aa_pattern = "aa"+code+task_id+"(03|09)\_05\.075\_v1\_3$"
    patterns = [aa_pattern]
    for suf in [".fasta", ".psipred", ".psipred_ss2"]:
        patterns.append(code+task_id+suf)
    
    try:
        for pattern in patterns:
            print "Pattern ",pattern
            
            # Find the results file in scr_job_dir
            match = paths.find(pattern, scr_job_dir)
            
            # If a results file matching pattern is not found in the scr_jobs_dir, raise exception.
            if not match:
                raise "Missing file"
            
            # Copy all results files in scr_job_dir to results_dir. 
            for path in match:
                # If the results file is a frag library, format the frag library file.
                if pattern == aa_pattern:
                    print "Filtering Columns of file %s" % path
                    # fileinput.input opens a file for iteration. Defining inplace=1 allows to given file to be
                    # altered "in place": file is backed up, and stdout is then directed to the file. 
                    # Filters the file: writes " " + first 47 chars of the line. (Specific formatting).
                    for line in fileinput.input(path, inplace=1):
                        print " "+line[:47].strip()
                    fileinput.close()
                print "Path ", path
                
                # Copy the result file to results_dir.
                utilities.copy(path, results_dir)
                
                # If the results file is a frag library, bzip it in the results_dir.
                if pattern == aa_pattern:
                    dest = paths.join(results_dir, paths.getFile(path))
                    utilities.system("bzip2 %s"  % dest)
                    
    # If gathering the results files fails, delete the results_dir and exit.
    except:
        paths.removerf(results_dir)
        raise


# Removes temporary fasta input and results file and temporary directories. 
def cleanup():
    paths.removerf(scr_job_dir)

# A wrapper for replicate_fasta, to copy fasta input file into a scratch location.
def replicate():

    # Create scr_job_dir
    paths.ensure(scr_job_dir)
    replicate_fasta()

# Copies input fasta file to the scr_fasta scratch directory (scr_job_dir/fasta_filename)
def replicate_fasta():
    
    if not paths.exists(scr_fasta):
        utilities.debug("Replicating fasta file")
        utilities.copy(fasta_file, scr_fasta)
        
    paths.existsOrFail(scr_fasta)


## Main functionality.

def main():
    global code, task_id, fasta_file, results_dir, scr_dir, scr_job_dir_dir, home, tools_dir
    
    # Create a parser and add command-line options.
    parser = OptionParser(usage=usage_str)
    set_options(parser)
    
    # Parse command-line options and arguments (does automatically on sys.argv[1:]).
    (options, args) = parser.parse_args()
    
    # Parse required positional arguments (stored in args via parse_args()).
    parse_positional_args(args, parser)

    # If task_id is specified, set to a 3-digit integer. If not, try to grab task_id from PBS env. var.
    if (options.task_id):
        task_id = '%(#)03d' % {"#": options.task_id}
    else:
        try:
            task_id = '%(#)03d' % {"#": int(os.environ['PBS_ARRAYID'])-1}
        except:
            raise Exception("No Task ID specified or in system.")
    
    # If properties file is specified, print usage and exit (not supported).
    if (options.properties_file):
        print "Specifying a properties file not currently supported. See usage."
        parser.print_help()
        sys.exit(1)
    
    
    # REPLACE THIS with getting fasta file to work on from tasks.pickle. How?
    # Replace '%s' in input fasta_file with code+task_id.
    try:
        fasta_file = fasta_file.replace("%s", code+task_id)
    except:
        print "Error: fasta_file {0} does not contain required \%s.".format(fasta_file)
        raise
    
    # Set results_dir to format: <results_dir>/<code>/<code>+<task_id>
    results_dir = os.path.join(results_dir, code, code+task_id)

    # Set job scratch directory to format: <scr_dir>/job/<code>+<task_id>
    scr_job_dir = os.path.join(scr_dir, "job", code+task_id)

    # TODO: Get rid of all this scr_fasta stuff, and just do the copy (see notes).
    # Set scr_fasta to <scr_job_dir><fasta_file filename only>
    scr_fasta = os.path.join(scr_job_dir, paths.getFile(fasta_file))

    # Sets home to CWD and tools_dir to CWD/tools. Hardcoded.
    # TODO: when make_fragments.local is incorporated, no tools will be used this way (instead, they'll all come from the Env.)
    home = os.getcwd()
    tools_dir = paths.join(home, "tools")
    
    # Run main fragmenting functionality.
    try:
        # Copy fasta input to scratch dir
        replicate()
        
        # Run fragment making script
        runScript(no_homs = options.nohoms)
        
        # Copy results to results dir
        gather()
        
    finally:
        if (options.cleanup):
            cleanup()
    
    utilities.debug("FINISHED!")
    
    
if __name__ == "__main__":
    main()


