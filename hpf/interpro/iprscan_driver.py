#!/usr/bin/env python

# A script to wrap executing IPRSCAN on many input files (a split fasta, for eg.)

import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description="A script to wrap IPRSCAN execution over many input files")
parser.add_argument("-d", "--directory", action="store", dest="in_dir", required=True,
        help="The directory containing fasta file inputs for IPRSCAN. Will use all files in dir")
parser.add_argument("-o", "--out_dir", action="store", dest="out_dir", required=True,
        help="Directory to put results files in")
parser.add_argument("-e", "--executable", action="store", dest="iprscan", default='iprscan',
        help="Executable path to IPRSCAN. Default is 'iprscan' (ie, iprscan is on PATH)")
parser.add_argument("-p", "--pfam_only", action="store_true", dest="pfam_only", default=False,
        help="Flag to execute iprscan with only PFAM method")
args = parser.parse_args()


in_files = os.listdir(args.in_dir)
os.chdir(args.out_dir)
for infile in in_files:
    print "Running iprscan on fasta input: {0}".format(infile)
    outfile = infile + ".interpro.out"
    cmd = [args.iprscan, "-cli", "-i", args.in_dir+infile, "-o", outfile, "-format", "raw", "-seqtype", "p", "-verbose", "-iprlookup", "-goterms"]
    print "Command:\t{0}".format(cmd)
    if args.pfam_only:
        cmd.append("-appl")
        cmd.append("hmmpfam")
    subprocess.check_call(cmd)
    print "Running iprscan on {0} Complete. Results in {1}".format(infile, outfile)


