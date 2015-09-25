#!/usr/bin/env python

import os
import sys
import getopt
from hpf.runtime import runtime, Flag
from hpf.utilities import consume
from hpf.processing import processor, SYNCHRONOUS
import logging

def remcm(sequence_key):
    import sys
    import os
    from hpf.hddb import tunnel, kill_tunnel 
    tunnel(5)
    from hpf.mcm.process import McmReRun
    from hpf.mcm import McmDB
    with open("/scratch/kd35/db175/data.mammothDb") as handle:
        mcmdb = McmDB.load(handle)
    mcmdb.scop='1.75'
    mcmdb.list_file = "/scratch/kd35/db175/list.mammoth"
    # Add int keys to the mcmdb._dict (this is rough)
    for key in mcmdb._dict.keys():
        mcmdb._dict[int(key)] = mcmdb._dict[key]
    jobdir = "/scratch/kd35/tmp/mcm/jobs/%i" % int(sequence_key)
    try:
    # Create an McmReRun task and run it
        if not os.path.exists(jobdir):
            os.mkdir(jobdir)
        os.chdir(jobdir)
        log_file = "/scratch/kd35/tmp/mcm/jobs/%i/outlog" % int(sequence_key)
        out_handle = open(log_file,"w")
        runtime().set_log(stdout=out_handle, stderr=out_handle)
        logging.basicConfig(filename=log_file, filemode='a')
        runtime().debug(os.environ['HOSTNAME'])
        task = McmReRun(jobdir,mcmdb,int(sequence_key))
        task.run()
    finally:
    # Clean up: gzip results, save log file, remove all other tmp files
        import subprocess
        runtime().debug("finally")
        for file in os.listdir(jobdir):
        if "data.mammoth" in file: 
            #save and gzip
            subprocess.Popen([r"gzip","-f",file])
        elif "outlog" in file:
            runtime().debug("saving outlog")
        else:
            os.remove(os.path.join(jobdir,file))
    #os.rmdir(jobdir)
    kill_tunnel()

    
def main(*args):
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS))
    runtime().debug("Using processor",pool)
    pool.make_tasks(None)
    consume(pool.run(remcm))
    
def _do(argv):
    r = runtime()
    r.description("""
    remcm.py [-options]
    Re-MCM a set of sequence keys.  Expects a 'tasks.pickle'.
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing"))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])

