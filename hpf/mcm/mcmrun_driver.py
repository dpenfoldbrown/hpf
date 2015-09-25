#
# A tool to encapsulate and organize the run functionality of MCM, corresponding
# to the HPF DB system to keep track of MCM runs and avoid redundancy.
# 
# @auth: dpb 11/21/2012
#
# NOTE: does a locking fetch of what code to work on from the HPF DB (table mcmRun).
# This is cool, as it avoids the wild task pool and PBS/SGE/etc array job and env
# var setup we were previously using. Should set ginzu up this way.
#

import os
import re
import sys
import argparse
import traceback

# CONFIG
ROSETTA_PATHS = "/scratch/kd35/mcm/env/paths.txt"
MAMMOTH_LIST  = "/scratch/kd35/mcm/env/list.mammoth"
MAMMOTH_DATA  = "/scratch/kd35/mcm/env/data.mammothDb"
NR_DATABASE   = "/scratch/kd35/ncbi_db/env_nr"
GINZU_VERSION = 4
IGNORE_GINZU_VERSION = True

DBSTORE = True      # cmdline parameter?
CLEANUP = True      # cmdline parameter?
DEBUG   = True      # cmdline parameter?


def main(letter_code, version, host, work_base, results_dir, retry):
    """
    Organizes the main functionality of getting a McmRun/Foldable to work on, extracting and checking
    the rosetta decoy file, and running MCM on the decoys. Everything after fetching the mcmRun DB
    record is in a try-catch block, to handle setting DB record status as well as possible.
    """
    from datetime import datetime
    from hpf.mcm.new_mcm import MCM
    from hpf.hddb.db import Session, McmRun
    session = Session()

    # Get MCM/foldable record (whole code, eg ox971) to work on
    print "DRIVER: Fetching MCM Run record from DB..."
    mcm_record = fetch_mcmRun_record(session, letter_code, version, host, retry)
    if not mcm_record:
        raise Exception("No McmRun record of code {0}, version {1} available to run (batch may be complete!)".format(letter_code, version))
    print "DRIVER: Obtained MCM Run record {1} for foldable: {0}".format(mcm_record.prediction_code, mcm_record.id)
    
    # Set up environment: make working directory and change there
    try:
        work_dir = os.path.join(work_base, mcm_record.prediction_code)
        if not os.path.isdir(work_dir): os.mkdir(work_dir)
        os.chdir(work_dir)
        print "DRIVER: Working directory holding all temp files: {0}".format(work_dir)
    except Exception as e:
        set_error(session, mcm_record, e, error_str="Error creating environment")
        raise e

    # Find decoy file and unzip to working directory (Exception if can't be found in results dir)
    try:
        results_file = fetch_results_file(results_dir, mcm_record.prediction_code)
        if not results_file: 
            raise Exception("No results file for code {0} found in {1}".format(mcm_record.prediction_code, results_dir))
        decoy_file = "{0}.decoys".format(mcm_record.prediction_code)
        decompress_results(results_file, decoy_file)
        print "DRIVER: Found results file '{0}'. Unzipped to working directory".format(results_file)
    except Exception as e:
        set_error(session, mcm_record, e, error_str="Error processing results file")
        raise e
    
    # Initialize MCM object and set McmRun record start time
    print "DRIVER: Creating MCM object on denovo results file: '{0}'".format(decoy_file)
    try:
        mcm = MCM(code=mcm_record.prediction_code,
                  decoy_file=decoy_file,
                  work_dir=work_base,
                  rosetta_pathsfile=ROSETTA_PATHS,
                  mammoth_listfile=MAMMOTH_LIST,
                  mammoth_datafile=MAMMOTH_DATA,
                  ginzu_version=GINZU_VERSION,
                  ignore_ginzu_version=IGNORE_GINZU_VERSION,
                  nr_db=NR_DATABASE,
                  dbstore=DBSTORE,
                  cleanup=CLEANUP,
                  debug=DEBUG,
                 )
        mcm_record.started = datetime.now()
    except Exception as e:
        set_error(session, mcm_record, e, error_str="Error initializing MCM")
        raise e

    # Run MCM object
    print "DRIVER: Running MCM process (start time set)"
    try:
        mcm.run()
    except Exception as e:
        set_error(session, mcm_record, e, error_str="Error running MCM")
        traceback.print_exc(file=sys.stderr)
        raise e
    
    # On successful finish, set finished time and status (only complete if results in DB) and exit
    mcm_record.finished = datetime.now()
    if DBSTORE:
        mcm_record.status  = 'complete'
        mcm_record.comment = "MCM run successful"
    else: 
        mcm_record.status  = 'unprocessed'
        mcm_record.comment = "MCM run successful, results not saved to DB"
    session.flush()

    print "DRIVER: MCM run on code '{0}', McmRun ID {1} complete".format(mcm_record.prediction_code, mcm_record.id)


def set_error(session, record, exception, error_str="Error"):
    """Sets the error in the McmRun DB record and prints error messages"""
    record.status  = 'error'
    record.comment = "{0}: {1}".format(error_str, exception)
    session.flush()
    print "DRIVER: MCM run on code '{0}', McmRun ID {1} failed".format(record.prediction_code, record.id)
    print "DRIVER: Setting McmRun status to error, halting with exception:\n\n"

def decompress_results(results_file, target_file):
    """Uses gzip module to decompress result file into target file"""
    import gzip
    results_handle = gzip.open(results_file, 'rb')
    target_handle = open(target_file, 'w')
    target_handle.write(results_handle.read())
    target_handle.close()
    results_handle.close()

def fetch_results_file(directory, prediction_code):
    """
    Searches a listing of given directory for a filename matching (by re)
    the given code. Returns first matching filename, or None if no files match.
    """
    if not os.path.isdir(directory):
        raise Exception("Directory {0} invalid for finding results file".format(directory))
    files = os.listdir(directory)
    for f in files:
        if re.search("{0}".format(prediction_code), f):
            return os.path.join(directory, f)
    return None

def fetch_mcmRun_record(session, letter_code, version, host, redo_error):
    """
    First attempts to lock the hpf.mcmRun table in order to ensure this process gets
    a prediction code that no other mcmrun_driver process is working on simultaneously.
    
    When a lock is acheived (it's a blocking call), queries the mcmRun table for records
    of the letter code and version specified with a status of 'unprocessed' (if 
    redo_error is True, then also with status 'error').

    If no record, returns None. Otherwise, sets the record's status to 'running',
    sets the host field, unlocks the table, and returns the record.
    """
    from sqlalchemy import or_
    from hpf.hddb.db import McmRun
    from sqlalchemy.sql.expression import func
    session.execute("LOCK TABLES {0} WRITE".format(McmRun.__tablename__))
    if redo_error:
        ready_record = session.query(McmRun).filter(McmRun.prediction_code.like("{0}%".format(letter_code)))\
                                            .filter(McmRun.version==version)\
                                            .filter(or_(McmRun.status=='unprocessed', McmRun.status=='error'))\
                                            .order_by(func.rand())\
                                            .first()
    else:
        ready_record = session.query(McmRun).filter(McmRun.prediction_code.like("{0}%".format(letter_code)))\
                                            .filter(McmRun.version==version)\
                                            .filter(McmRun.status=='unprocessed')\
                                            .first()
    if ready_record:
        ready_record.status = 'running'
        ready_record.host = host
        session.flush()
    session.execute("UNLOCK TABLES")
    return ready_record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Driver tool to get McmRun record from DB (see mcmrun_setup) and run MCM on it")
    parser.add_argument("--code", action="store", dest="letter_code", required=True,
            help="The letter prediction code of the batch of foldable records to run MCM on. EG: oa")
    parser.add_argument("--version", action="store", type=int, dest="version", default=1,
            help="The version of the MCM run. Simply an arbitrary integer to allow for keeping track of multiple MCM runs, if desired")
    parser.add_argument("--host", action="store", dest="host", required=True,
            help="A string to identify where the driver is being run. Machine name, cluster name, etc.")
    parser.add_argument("--work_dir", action="store", dest="work_dir", default="/tmp/mcm",
            help="The base directory to create a working MCM directory (workdir/<code>) in")
    parser.add_argument("--results_dir", action="store", dest="results_dir", required=True,
            help="The directory containing the gzipped rosetta prediction files, eg: ox999.result.gz")
    parser.add_argument("--redo_error", action="store_true", dest="redo_error", required=False,
            help="If present, will cause driver to run MCM on McmRun DB entries with 'error' status")

    args = parser.parse_args()

    # Check given code
    if not re.match(r"[a-z]{2}", args.letter_code):
        raise Exception("Code {0} not in valid form (two lower-case letters)".format(args.letter_code))
    
    # Set and check directories
    args.work_dir = os.path.abspath(os.path.expanduser(args.work_dir))
    args.results_dir = os.path.abspath(os.path.expanduser(args.results_dir))
    if not os.path.isdir(args.work_dir):
        raise Exception("Working directory {0} is not valid".format(args.work_dir))
    if not os.path.isdir(args.results_dir):
        raise Exception("Results directory {0} is not valid".format(args.results_dir))

    main(args.letter_code, args.version, args.host, args.work_dir, args.results_dir, args.redo_error)

