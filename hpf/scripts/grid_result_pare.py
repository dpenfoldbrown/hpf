#!/usr/bin/env python

# A script to read in all directories representing denovo results from the grid,
# and for each file in those directories, gunzip the file, pare the results (keep
# top 25%), then rezip the file with bzip2 -9 compression, and delete the old file.

# To know which directories and files to pare, will parse a dictionary of
# letter-code -> # of codes not MCMed in letter.
# If # not MCMed is <300, all files in the letter-code dir will be pared
# If # not MCMed is >=300, only files that have been MCMed in that letter-code
# dir will be pared.
# (Requires two files).

import os
import re
import subprocess
from hpf.hddb.mcm.denovo_results import denovoResultFile


# Best percentage of pared files to keep
PERCENT_TO_KEEP = 25

# The directory containing grid results directories
#base_dir = '/scratch/dpb3/denovo_results'
base_dir = '/archive/rb133/hpf2/results'

# The hpf1 and hpf2 directory name formats for parsing letter and full prediction codes
hpf2_letcode_pat = r"(?P<letcode>[a-z]{2})_hpf2_results"
hpf1_letcode_pat = r"(?P<letcode>[a-z]{2})"

# List files for codes w/o MCM and counts for letter codes
nomcm_file = '/home/dpb3/pare/nomcm_codes.txt'
nomcm_count_file = '/home/dpb3/pare/nomcm_letcode_count.txt'

nomcm_codes = list()
handle = open(nomcm_file)
for line in handle:
    nomcm_codes.append(line.rstrip())
handle.close()
print "Number of codes without MCM results: {0}".format(len(nomcm_codes))

nomcm_count = dict()
handle = open(nomcm_count_file)
for line in handle:
    code, count = line.rstrip().split()
    nomcm_count[code] = int(count)
handle.close()
print "Compiling count dictionary complete"


def main():
    results_dirs = os.listdir(base_dir)
    for dir in results_dirs:
        try:
            letter_code = re.match(hpf2_letcode_pat, dir).group('letcode')
        except AttributeError as e:
            print "Letter code pattern for matching directory '{0}' failed".format(dir)
            raise(e)
        
        print "Processing {0} results".format(letter_code)
        
        results_dir = os.path.join(base_dir, dir)
        results_files = remove_nonresult_files( os.listdir(results_dir) )
        if letter_code not in nomcm_count.keys():
            print "Prediction code {0} not in set of counted prediction codes. Skipping..".format(letter_code)
            continue
        elif nomcm_count[letter_code] >= 300:
            os.chdir(results_dir)
            pare_some(results_files)
        else:
            os.chdir(results_dir)
            pare_all(results_files)
    print "Reducing grid results storage complete"


def test():
    results_dirs = os.listdir(base_dir)
    for dir in results_dirs:
        try:
            letter_code = re.match(hpf2_letcode_pat, dir).group('letcode')
        except AttributeError as e:
            print "Letter code pattern for matching directory '{0}' failed".format(dir)
            raise(e)
        
        print "Testing {0} results".format(letter_code)
        
        results_dir = os.path.join(base_dir, dir)
        results_files = remove_nonresult_files( os.listdir(results_dir) )
        os.chdir(results_dir)
        for file in results_files:
            test_parse(file)
    print "Testing all files complete"

def remove_nonresult_files(files):
# Takes a list of files, returning only those matching the filename format of un-pared denovo results files
    hpf2_filepat = r"[a-z]{2}[0-9]{3}\.result\.gz"
    hpf1_filepat = r"[a-z]{2}[0-9]{3}\.out\.gz"
    results_files = list()
    for file in files:
        if re.match(hpf2_filepat, file) or re.match(hpf1_filepat, file):
            results_files.append(file)
    return results_files

def get_prediction_code(filename, hpf1=None, hpf2=None):
# Returns a prediction code string from a filename string
    if hpf1:
        hpf1_predcode_pat = r"(?P<predcode>[a-z]{2}[0-9]{3})\.out\.gz"
        try:
            prediction_code = re.match(hpf1_predcode_pat, filename).group('predcode')
        except AttributeError as e:
            print "Prediction code pattern for matching file '{0}' failed".format(file)
            raise(e)
    else:
        hpf2_predcode_pat = r"(?P<predcode>[a-z]{2}[0-9]{3})\.result\.gz"
        try:
            prediction_code = re.match(hpf2_predcode_pat, filename).group('predcode')
        except AttributeError as e:
            print "Prediction code pattern for matching file '{0}' failed".format(file)
            raise(e)
    return prediction_code

def pare_some(results_files):
# Pare only those results that are NOT in list of codes that have not been MCMed
    
    print "(pare_some)"

    for file in results_files:
        prediction_code = get_prediction_code(file, hpf2=True)
        if prediction_code in nomcm_codes:
            print "Code {0} has not been MCMed. Skipping..".format(prediction_code)
            continue
        pare(file, prediction_code=prediction_code)

def pare_all(results_files):
# Pare all results in given list
    
    print "(pare_all)"
    
    for file in results_files:
        pare(file)

def pare(file, prediction_code=None):
# Gunzip, parse, cut down to 25% best, bzip2 -9 pared file, delete old file
    if not prediction_code:
        prediction_code = get_prediction_code(file, hpf2=True)
    
    print "\tParing results for {0}".format(prediction_code)

    # Gunzip given file (cut '.gz' off the filename)
    if subprocess.call(["gunzip", file]) != 0:
        raise Exception("Gunzipping file '{0}' failed".format(file))
    decompressed_file = file.rstrip('.gz')

    # Open, parse, and pare file, storing in new file. Via denovo_results.denovoResultFile obj
    drf = denovoResultFile(decompressed_file, prediction_code)
    pared_file = decompressed_file + ".reduced"
    drf.write_to_file(outfile=pared_file, percent=PERCENT_TO_KEEP)

    # Bzip2 -9 new pared file
    if subprocess.call(["bzip2", "-9", pared_file]) != 0:
        raise Exception("Bzip2 -9'ing file '{0}' failed".format(pared_file))
    recompressed_file = pared_file + ".bz2"

    # Delete old file (NOTE: DO NOT DELETE FOR NOW)
    if subprocess.call(["rm", decompressed_file]) != 0:
        raise Exception("Deleting old file '{0}' failed".format(decompressed_file))
    #print "Would delete old file '{0}' here..".format(decompressed_file)
    
    print "\tParing '{0}' complete".format(file)

def test_parse(file, prediction_code=None):
# Test the denovo result file parsing without doing any other processing
    if not prediction_code:
        prediction_code = get_prediction_code(file, hpf2=True)
    
    print "\tTesting results for {0}".format(prediction_code)

    # Gunzip given file (cut '.gz' off the filename)
    if subprocess.call(["gunzip", file]) != 0:
        raise Exception("Gunzipping file '{0}' failed".format(file))
    decompressed_file = file.rstrip('.gz')

    # Open, parse, and pare file, storing in new file. Via denovo_results.denovoResultFile obj
    drf = denovoResultFile(decompressed_file, prediction_code)

    # Re-Gzip file
    if subprocess.call(["gzip", decompressed_file]) != 0:
        raise Exception("Re-Gzipping file '{0}' failed".format(decompressed_file))

    print "\tTesting result {0} completed".format(prediction_code)


if __name__ == "__main__":
    #main()
    test()

