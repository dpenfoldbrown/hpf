#!/usr/bin/env python

# Script to drive a sequential fragmentor run (via hpf.fragmentor.fragmentor). Could
# attempt to use task pool system, but for frag picking, this results in race condi-
# tions for creating output files, causing success checks to fail. So we'll avoid
# that.
# dpb, 8/10/2012

import os
import re
import argparse
from glob import glob
from hpf.fragmentor.fragmentor import Fragmentor


def frag_driver(code, fasta_files, results_dir, config_file, nohoms, cleanup):
    """Will loop through the given fasta files, creating and running a Fragmentor obj for each"""
    print "sequential frag:: Running Sequential Fragmentor on {0} files".format(len(fasta_files))
    success = 0
    for file in fasta_files:
        print "sequential frag:: Running fragmentor on {0}".format(file)
        f = Fragmentor(fasta_file=file, code=code, results_dir=results_dir, config_file=config_file)
        try:
            f.run(nohoms=nohoms, cleanup=cleanup)
            success += 1
        except Exception as e:
            print "sequential frag:: Fragmentor run on {0} failed".format(file)
            raise e
    print "sequential frag:: Fragmentor run successfully on {0} fastas".format(success)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("A driver for sequentially running Fragmentor on fasta files")
    parser.add_argument("-n", "--nohoms", action="store_true", dest="nohoms", default=False,
            help="Flag. Specifies passing -nohoms to the fragment picker (no homologues)")
    parser.add_argument("-c", "--no_cleanup", action="store_false", dest="cleanup", default=True,
            help="Flag. add -c to NOT cleanup working files")
    parser.add_argument("-f", "--config_file", action="store", dest="config_file", default=None,
            help="The config file to pass to the Fragmentor object. Default None (looks in hpf.fragmentor dir)")
    parser.add_argument("--code", action="store", dest="code", required=True, 
            help="The two-letter code of the fasta files to run")
    parser.add_argument("--input", action="store", dest="fasta_in", required=True,
            help="Input file string representing fasta files. Can be a single file or a RE expanding to multiples files")
    parser.add_argument("--results_dir", action="store", dest="results_dir", required=True,
            help="The directory to store results in")
    args = parser.parse_args()

    # Check parameters
    if not re.match(r"[a-z]{2}", args.code):
        raise Exception("Code {0} not valid".format(args.code))
    
    if not os.path.isdir(args.results_dir):
        raise EnvironmentError("Results directory {0} not valid".format(args.results_dir))

    if args.config_file and not os.path.isfile(args.config_file):
        raise EnvironmentError("Config file {0} not valid".format(args.config_file))
    
    fasta_files = glob(os.path.expanduser(args.fasta_in))
    for f in fasta_files:
        if not os.path.isfile(f):
            raise EnvironmentError("Expanded file {0} not valid (use full paths)".format(f))


    frag_driver(args.code, fasta_files, args.results_dir, args.config_file, args.nohoms, args.cleanup)
