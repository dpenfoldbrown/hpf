"""
Extractor tools (used in MCM) based on the Rosetta++ rosetta.gcc -extract tool
Auth.: dpb 4/05/2012
"""

import os
import subprocess

# Example execution:
# pdbfilename = extract(decoy_file, 'S_0027_1715', log_file=log, paths=pathfile, debug=True)

def extract(source_file, tag, log_file=None, paths=None, bin="rosetta", debug=False):
    """A function to extract the PDB atom record of a decoy in silent form in source_file
    described by the Rosetta "description" (S_id) 'tag'.
    Params:
        source_file - the silent file of decoys to extract 'tag' decoy from
        tag   - a Rosetta "description" of the form S_1234_1234 in the given source_file
        log_file    - a file to dump STDOUT and STDERR into
        paths - the location of the Rosetta paths.txt file, required in CWD or explicitly given
        bin   - the binary to execute
    Returns filename of extracted decoy record
    """
    if debug: print "Extract::"
    if not os.path.isfile(source_file):
        raise Exception("Given source file '{0}' invalid".format(source_file))

    # Build extractor command
    cmd = [bin,]
    if paths: cmd += ["-paths", paths]
    cmd += ["-extract", "-s", source_file, "-t", tag]
    if debug: print "Command: ", cmd

    # Run extractor command
    if log_file:
        with open(log_file, 'a') as log_handle:
            subprocess.check_call(cmd, stdout=log_handle, stderr=subprocess.STDOUT)
    elif debug:
        subprocess.check_call(cmd)
    else:
        subprocess.check_call(cmd, stdout=os.devnull, stderr=subprocess.STDOUT)

    # Rosetta extractor will output pdb in this format. Check for file
    pdb_file = "{0}.pdb".format(tag)
    if not os.path.isfile(pdb_file):
        raise OSError("Expected extractor output file '{0}' missing. Extract failed".format(pdb_file))
    if debug: print "Extractor produced pdb file '{0}'".format(pdb_file)
    
    return pdb_file

def extract_all():
    assert 0, "Not yet implemented"

def extract_list():
    assert 0, "Not yet implemented"

