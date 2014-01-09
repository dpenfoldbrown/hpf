#!/usr/bin/env python

# A simple script to split a multiple-sequence fasta file into a number of different files.
#
# - Open fasta file and parse sequences into record list
# - Calc number of sequences per split file
# - Loop through sequences, adding to file. Create new file every mod result = 0

import argparse
from os import chdir
from os.path import isdir
from Bio import SeqIO

# Commandline arguments:
parser = argparse.ArgumentParser(description="A fasta file splitter")
parser.add_argument("-f", "--fasta", action="store", dest="fasta_file", required=True, 
        help="The fasta file to split")
parser.add_argument("-o", "--out_prefix", action="store", dest="outfile_prefix", required=True, 
        help="The prefix of the outfiles")
parser.add_argument("-d", "--out_dir", action="store", dest="out_dir", default=".", 
        help="The directory to output files to (default is current dir)")
parser.add_argument("-p", "--num_pieces", action="store", dest="num_pieces", required=True,
        help="The number of pieces to split fasta file into")
args = parser.parse_args()

# Change working directory (where files will be output)
chdir(args.out_dir)

# Open file and parse with BioPython
handle = open(args.fasta_file, "rU")
seqs = list(SeqIO.parse(handle, "fasta"))
handle.close()

num_seqs = len(seqs)
seqs_per_file = num_seqs / int(args.num_pieces)
file_num = 0

for i in range(num_seqs):
    if (i % seqs_per_file == 0):
        if i != 0:
            handle.close()
        filename = "{0}{1}".format(args.outfile_prefix, file_num)
        print "Creating new file '{0}' containing sequences {1} to {2} of {3}".format(filename, file_num*seqs_per_file, file_num*seqs_per_file+seqs_per_file-1, num_seqs)
        handle = open(filename, "w")
        file_num += 1

    #print ">{0}\n{1}".format(seqs[i].description, seqs[i].seq)
    handle.write(">{0}\n{1}\n".format(seqs[i].description, seqs[i].seq))

handle.close()    
print "Splitting fasta complete"

