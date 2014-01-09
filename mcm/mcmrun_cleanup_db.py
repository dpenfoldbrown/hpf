#
# A sanity check tool to cleanup mcm runs which may not have finished properly 
# according to the HPF DB system
# 
# @auth: kdrew 2/13/2013
#
# NOTE: checks mcm tables for completeness and updates HPF DB (table mcmRun).
#	revised from mcmrun_driver.py
# 
# 
#

import os
import re
import argparse

# CONFIG

def main(letter_code, version):
    """
    Organizes the main functionality of getting a McmRun/Foldable to work on, extracting and checking
    the rosetta decoy file, and running MCM on the decoys. Everything after fetching the mcmRun DB
    record is in a try-catch block, to handle setting DB record status as well as possible.
    """
    from datetime import datetime
    from hpf.mcm.new_mcm import MCM
    from hpf.hddb.db import Session, McmRun, RosettaConvergence, RosettaCluster, Structure, Mammoth, McmData
    session = Session()

    # Get MCM/foldable record (whole code, eg ox971) to work on
    print "DRIVER: Fetching MCM Run record from DB..."
    mcm_records = fetch_mcmRun_records(session, letter_code, version)
    if not mcm_records:
        raise Exception("No McmRun records of code {0}, version {1} available".format(letter_code, version))
    for mcm_record in mcm_records:
        print "CLEANUP: Obtained MCM Run record {1} for foldable: {0}, status {2}".format(mcm_record.prediction_code, mcm_record.id, mcm_record.status)

        #kdrew: check for proper values in mcm tables
        rosetta_convergence_count = 0 #1 entry per mcm
        rosetta_cluster_count = 0 #normally 25 entries per mcm
        structure_count = 0 #normally 25 entries per mcm, should equal rosetta_cluster_count
        mammoth_count = 0 #normally 5 entries
        mcm_count = 0 #normally 5 entries, should equal mammoth_count

        rosetta_convergence_error = False
        rosetta_cluster_error = False
        structure_error = False
        mammoth_error = False
        mcm_error = False

        error_str = ""

        rosetta_convergence_count = session.query(RosettaConvergence).filter_by(outfile_key = mcm_record.foldable.id).count()
        rosetta_convergence = session.query(RosettaConvergence).filter_by(outfile_key = mcm_record.foldable.id).first()
        print rosetta_convergence, rosetta_convergence_count

        if rosetta_convergence:
            rosetta_cluster_count = session.query(RosettaCluster).filter_by(convergence_key=rosetta_convergence.id).count()
            rosetta_clusters = session.query(RosettaCluster).filter_by(convergence_key=rosetta_convergence.id).all()
            print rosetta_clusters[0], len(rosetta_clusters)
            if len(rosetta_clusters) == 0:
                rosetta_cluster_error = True
                error_str = error_str + " missing rosetta_cluster entry,"

            for rosetta_cluster in rosetta_clusters:
                structure = session.query(Structure).filter_by(id=rosetta_cluster.structure_key).first()
                #print structure
                if structure:
                    structure_count += 1
                else:
                    structure_error = True
                    error_str = error_str + " missing structure entry,"
        else:
            rosetta_convergence_error = True
            error_str = error_str + " missing rosetta_convergence entry,"
               
        print "structure_count: %s" % (structure_count,)

        mcm_entries = session.query(McmData).filter_by(outfile_key = mcm_record.foldable.id).all()
        if mcm_entries:
            mcm_count = len(mcm_entries)
            for mcm_entry in mcm_entries:
                print mcm_entry
                mammoth_entry = session.query(Mammoth).filter_by(p_structure_key = mcm_entry.structure_key).all()
                if mammoth_entry:
                    mammoth_count += 1
                else:
                    mammoth_error = True
                    error_str = error_str + " missing mammoth entry,"
        else:
            mcm_error = True
            error_str = error_str + " missing mcm entry,"

        print "mammoth_count: %s" % (mammoth_count,)

        #kdrew: if entries are in the db and no errors
        if not mammoth_error and not mcm_error and not rosetta_convergence_error and not rosetta_cluster_error and not structure_error:
            #kdrew: if counts are as expected
            if mammoth_count == 5 and mcm_count == 5 and rosetta_cluster_count == 25 and structure_count == 25:
                if mcm_record.status != 'complete': 
                    mcm_record.finished = datetime.now()
                    mcm_record.status  = 'complete'
                    mcm_record.comment = "MCM run successful (on mcmrun_cleanup_db)"
                    session.flush()
                #else: do not change completed runs
            #kdrew: missing db entries but no errors
            else:
                warning_str = ""
                if mammoth_count != 5:
                    warning_str = warning_str + " missing mammoth entries,"
                if mcm_count != 5:
                    warning_str = warning_str + " missing mcm entries,"
                if rosetta_cluster_count != 25:
                    warning_str = warning_str + " missing rosetta_cluster entries,"
                if structure_count != 25:
                    warning_str = warning_str + " missing structure entries,"

                mcm_record.finished = datetime.now()
                mcm_record.status  = 'complete'
                mcm_record.comment = "MCM run successful (on mcmrun_cleanup_db), number of entries in db not as expected %s" % (warning_str,)
                session.flush()
        #kdrew: error occurred
        else:
            if mcm_record.status != 'unprocessed':
                mcm_record.finished = datetime.now()
                mcm_record.status  = 'error'
                mcm_record.comment = "MCM run not successful (on mcmrun_cleanup_db), %s" % (error_str,)
                session.flush()



        print "CLEANUP: MCM run on code '{0}', McmRun ID {1} ".format(mcm_record.prediction_code, mcm_record.id)
    


def fetch_mcmRun_records(session, letter_code, version):
    """
        Fetches all mcmRun records with prediction code and version
    """
    from hpf.hddb.db import McmRun
    ready_records = session.query(McmRun).filter(McmRun.prediction_code.like("{0}%".format(letter_code)))\
                                        .filter(McmRun.version==version)\
                                        .all()
    return ready_records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sanity check tool to update McmRun record from DB (see mcmrun_setup) and check for completion")
    parser.add_argument("--code", action="store", dest="letter_code", required=True,
            help="The letter prediction code of the batch of foldable records to run MCM on. EG: oa")
    parser.add_argument("--version", action="store", type=int, dest="version", default=1,
            help="The version of the MCM run. Simply an arbitrary integer to allow for keeping track of multiple MCM runs, if desired")

    args = parser.parse_args()

    # Check given code
    if not re.match(r"[a-z]{2}", args.letter_code):
        raise Exception("Code {0} not in valid form (two lower-case letters)".format(args.letter_code))
    
    main(args.letter_code, args.version)

