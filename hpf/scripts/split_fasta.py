#!/usr/bin/env python

import sys
import getopt
import subprocess
import os
import math
from Bio.SeqIO.FastaIO import FastaWriter
from Bio import SeqIO
from hpf.runtime import runtime, IntegerOption, Flag
from hpf.processing import processor, MultiProcessor, SYNCHRONOUS
from hpf.utilities import consume

PARTS = "parts"

def split(fasta_file, parts):
    cmd = "grep -c '>' %s"%fasta_file
    out,err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
    count = int(out.strip())
    part_size = int(math.ceil(float(count)/float(parts)))
    runtime().debug("Part size",part_size)
    writer = None
    handle = None
    with open(fasta_file) as fasta:
        for i,record in enumerate(SeqIO.parse(fasta,"fasta")):
            if i%part_size==0:
                part = i/part_size
                newfile=fasta_file+".%i"%part
                runtime().debug(i,i%part_size,newfile)
                if handle:
                    handle.close()
                handle = open(newfile,"w")
                #print handle
                if writer:
                    writer.write_footer()
                writer =  FastaWriter(handle)
                writer.write_header()
            #print record
            writer.write_record(record)
    
def main(*args):
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS))
    runtime().debug("Using processor",pool)
    pool.make_tasks(lambda: [(fasta, runtime().opt(PARTS)) for fasta in args])
    consume(pool.run(lambda args: split(*args)))
    
def _do(argv):
    r = runtime()
    r.description("""
    split_fasta.py [-options] *fastafiles
    Splits each fasta into the given number of parts.
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing", default=True))
    r.add_option(IntegerOption(PARTS, "p", description="Number of parts to split the fasta into", default=20))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
