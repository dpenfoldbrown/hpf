#!/usr/bin/python

"""
A Quick script to drive importing HHPred results into the DB. Made quick and dirty, rewrite as you please
dpb 5/25/2012
"""

import os
import re
from hpf.hhpred import HpfHHPredResultFile


# Add all "matching" files in directories in list
directories = list()

directories.append('/Users/dpb/bonneau-dev/sandbox/hhpred')
directories.append('/Users/dpb/bonneau-dev/sandbox/hhpred2')

#directories.append('/scratch/kd35/hh_search_runs/human/search')
#directories.append('/scratch/kd35/hh_search_runs/humanRNA/search')
#directories.append('/scratch/kd35/hh_search_runs/mouse/search')

# Files must match this pattern to be parsed. (COULD make this capture DB and seqkey)
FILE_PATTERN = r"[0-9]+_scop175\.hhr"

# Database string and version (see hpf.hhpred_version DB table) to pass to object
DB = "scop_1.75"
VERSION = 1

# Import by directory
for dir in directories:
    print "DRIVER -- Importing HHPred results found in '{0}'".format(dir)
    imported = 0
    skipped = 0
    files = os.listdir(dir)
    for file in files:
        if re.match(FILE_PATTERN, file):
            print "DRIVER -- Importing results from file '{0}'".format(file)
            hhprf = HpfHHPredResultFile(results_file=os.path.join(dir, file), database=DB, version_key=VERSION)
            hhprf.db_store_hits()
            imported += 1
        else:
            skipped += 1
            continue
    print "DRIVER -- Imported {0} files, skipped {1} from directory '{2}'".format(imported, skipped, dir)
print "DRIVER -- Importing HHPred results for all directories complete"
