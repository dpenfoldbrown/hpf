#!/usr/bin/env python

## A python fragmentor for the creation of fragment libraries and other analyses from sequences. 
## Uses nnmake, psipred, and other alignment tools to support fragmenting functionality. (nnmake 
## ultimately responsible for fragment library creation).
##
## Essentially a rewrite and incorporation of make_fragments.local.pl and frag_array.py, themselves
## wrappers and drivers for fragment picking functionality.
##
## Auth: Duncan Penfold-Brown, 1/10/2011. 
## From frag_array.py (Patrick Winters) and make_fragments.local (rohl, bonneau, roberston 7/25/2002).
##
## Also serves as an example driver for parallelization

## Imports
import sys
from glob import glob
from os.path import expanduser, isdir

from hpf.fragmentor import fragmentor
from optparse import OptionParser
from hpf.utilities import consume
from hpf.processing import processor


## Globals

usage_str = """Usage: %prog [-n] [-c] [-f CONFIG FILE]  CODE  FASTA_FILE  RESULTS_DIR
  
  CODE          The two-letter prediction code (used in naming results files). 
                  Will automatically be chopped to first two given chars.
  FASTA_FILE    Input file string representing fasta files. Can be a single file 
                  or a RE expanding to multiples files (eg. "/dir/*.fasta"). If
                  using a RE for expanding a directory, PUT IN QUOTES (so the RE
                  (is not expanded automatically by the shell).
  RESULTS_DIR   The directory for storing fragmentor results. Results are moved 
                  to RESULTS_DIR/[CODE]/[CODE][#].
"""

code = None
fasta_in = None
results_dir = None
nohoms = None
cleanup = None
frag_config = None


## Option Setting and Parsing.

# Sets the command-line options for a passed-in OptionParser object (via optparse, Python 2.3+)
def set_options(parser):
    parser.add_option("-n", action="store_true", dest="nohoms", default=False,
                      help="Flag. Specifies -nohoms to the fragment picker.")
    parser.add_option("-c", action="store_false", dest="cleanup", default=True,
                      help="Flag. Do NOT clean up. Defaults to TRUE (by default, cleanup will occur).")
    parser.add_option("-f", dest="config_file", metavar="FILE",
                      help="The config file to pass to Fragmentor object, specifying tools and databases.")


# Parses positional arguments necessary for functionality (CODE, FASTA_INPUT, RESULTS_DIR)
# Takes a list of args and a parser Option object.
def parse_positional_args(args, parser):
    if len(args) < 3:
        print "Error: Incorrect usage."
        parser.print_help()
        sys.exit(1)
    
    # Return code, fasta_in, and results_dir to main().
    # TODO: Do some verifcation on positional inputs - check for file existence, etc.
    code = str(args[0])[:2]
    fasta_in = args[1]
    results_dir = args[2]
    #DEBUG
    print "parse_positional_args::\n\tcode: {0}\n\tfasta_in: {1}\n\tresults_dir: {2}".format(code, fasta_in, results_dir)
    
    if not isdir(results_dir):
        print "Error: RESULTS_DIR field must be a directory. Current value: {0}".format(results_dir)
        parser.print_help()
        sys.exit(1)
    
    return (code, fasta_in, results_dir)


## Fragmentor and Task functionality

## Runs the fragment picking functionality on the given fasta file/sequence:
#   Sets up and creates working directory structure; copy input fasta to working directory; run 
#   fragmenting functionality tool chain; filter final results and move to results dir; clean up
#   working directory.
# Called by processing_pool.run
def frag_driver(fasta_file):
    print "frag_driver:: Fasta file: {0}.".format(fasta_file)
    
    # Create a fragmentor object (via hpf.fragmentor) and call it's run function.
    frag = fragmentor.Fragmentor(fasta_file=fasta_file, code=code, results_dir=results_dir, config_file=frag_config)
    try:
        frag.run(nohoms=nohoms, cleanup=cleanup)
    except:
        print "Fragmentor run on {0} failed.".format(fasta_file)
        raise


## Creates and returns a list of tasks, where a task is a fasta input file derived from given input
#   string (fasta_in, of the form /somedir/*.fasta or /somedir/onefile.fasta). Must be a string that
#   can be glob'ed to return a file or number of files (eg, a regex representing multiple files).
# Called by  processing_pool.make_tasks.
def tasks(fasta_in):
    print "tasks:: Input string: {0}.".format(fasta_in)
    
    # Parse fasta input string (single file or file regexp) into a list of files rep'ing frag tasks.
    fasta_files = glob(expanduser(fasta_in))
    if fasta_files == []:
        raise Exception("Input string {0} contains no valid files.".format(fasta_in))
    print "tasks:: Input fasta files: ", fasta_files
    
    # Return list of (fasta) files representing tasks to send to the 'frag_driver' function.
    return fasta_files


## Main functionality: create a parser and read command-line options and positional params; 
#   set up a processing pool and serialize execution tasks; run tasks through frag_driver function
#   on processing pool.
def main():
    global code, fasta_in, results_dir, nohoms, cleanup, frag_config

    # Create a parser, add options, parse options, and parse leftover required positional args.
    parser = OptionParser(usage=usage_str)
    set_options(parser)
    (options, args) = parser.parse_args()
    (code, fasta_in, results_dir) = parse_positional_args(args, parser)
    
    # Set config file, if specified.
    if options.config_file:
        frag_config = options.config_file
    
    # Set remaining globals
    nohoms = options.nohoms
    cleanup = options.cleanup
    
    # Create a processing pool (via hpf.processing) and make tasks to serve to the frag_driver.
    print "main:: Creating processor pool."
    pool = processor()
    print "main:: Serializing tasks based on input string: {0}.".format(fasta_in)
    pool.make_tasks(tasks, fasta_in)
    print "main:: Running tasks."
    consume(pool.run(frag_driver))


## Main Call.
if __name__ == "__main__":
    main()
