#!/usr/bin/env python

# A script to package fragmentor results for IBM/WCG
# Must be run in the directory where code directories exist
# (eg. if packaging .../or/, .../os, must be in ...)

import os
import sys
import subprocess

# MD5 exec. (changes from md5 to md5sum)
md5_exe = "md5sum"

# Manually populate directory and list of codes to package.
dir = "/scratch/dpb3/fragmentor/results"
codes = ["ql", "qm", "qn", "qo", "qp"]

os.chdir(dir)

for code in codes:
    # Check for results directory existence.
    if not os.path.exists(code):
        print "Error: directory {0} does not exist in the current directory".format(code)
        sys.exit(1)
    
    # Get results range.
    results_dirs = sorted(os.listdir(code))
    low  = results_dirs[0][2:]
    high = results_dirs[-1][2:]
    print "Packaging results in {0}, containing {3} results, {1} - {2}".format(code, low, high, len(results_dirs))
    
    # Create and execute a command to tar and bzip results dir.
    compress_name = "{0}.{1}-{2}.tar.bz2".format(code, low, high)
    cmd = ["tar", "cjf", compress_name, code]
    print "Compression command: ", cmd
    print "Compressing results directory..."
    if subprocess.call(cmd) != 0:
        print "Error: compression command failed: ", cmd
        sys.exit(1)

    # Take an md5 sum of the compressed results and store it in file.
    md5_name = "{0}.{1}-{2}.tar.md5".format(code, low, high)
    cmd = "{0} {1} > {2}".format(md5_exe, compress_name, md5_name)
    print "MD5 command: {0}".format(cmd)
    print "Computing MD5 sum of {0} file...".format(compress_name)
    if os.system(cmd) != 0:
        print "Error: md5 sum command failed: {0}".format(cmd)
        sys.exit(1)

    print "Packaging for result {0} is complete.".format(code)

print "Packaging results for following codes complete: ", codes
 
