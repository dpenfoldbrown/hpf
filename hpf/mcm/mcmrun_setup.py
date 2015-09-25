#
# A setup tool to initialize the HPF DB for an MCM run. 
# Adds and inits entries for each foldable prediction code to the hpf.mcmRun table,
# while making checks on data already in DB (to avoid redundancy)
#
# @auth: dpb 11/21/2012
#
# NOTE: version (here and in the HPF DB) could be overwritten to support a more
# robust and legit version number/system, such as that for the ginzuRun ginzu_version
#

import re
import argparse

def mcm_setup(letter_code, version, comment):
    """
    Takes a foldable record batch prediction code (two letters, eg 'ax'), a version 
    (which serves as an arbitrary identifier for keeping track of MCM runs), and a comment
    to add the to hpf.mcmRun row's comment field.
    
    Checks DB for code and version. If pre-existing, exits. Otherwise, creates mcmRun
    entries corresponding to each foldable code, to keep track of their MCM processing.
    See mcmrun_driver.py for the actual MCM fetch and run functionality.
    """
    from datetime import datetime
    from hpf.hddb.db import Session, McmRun, FilesystemOutfile as FoldableRecord, push_to_db
    session = Session()

    # Check code
    if not re.match(r"[a-z]{2}", letter_code):
        raise Exception("Code {0} does not match code format (two lower-case letters)".format(code))

    # Check for pre-existing McmRun records with same code and version
    if session.query(McmRun).filter(McmRun.prediction_code.like("{0}%".format(letter_code))).filter(McmRun.version==version).first():
        raise Exception("Code {0}, version {1} already exists in DB. Use a different version or clean DB").format(letter_code, version)

    # Get all foldables, and create McmRun objects from them and passed in params
    foldables = session.query(FoldableRecord).filter(FoldableRecord.prediction_code.like("{0}%".format(letter_code))).all()
    insert_date = datetime.now()
    mcmrun_count = 0
    for fold_rec in foldables:
        new_mcmrun = McmRun(prediction_code=fold_rec.prediction_code,
                            sequence_key=fold_rec.sequence_key,
                            inserted=insert_date,
                            version=version,
                            comment=comment,
                           )
        push_to_db(session, new_mcmrun, exception_str="Could not push McmRun code {0} to DB".format(fold_rec.prediction_code))
        mcmrun_count += 1

    print "McmRun setup: {0} McmRun objects pushed to DB".format(mcmrun_count)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Setup tool to initialize an MCM run in the HPF DB")
    parser.add_argument("--code", action="store", dest="letter_code", required=True,
            help="The letter prediction code of the batch of foldable records to init MCM on")
    parser.add_argument("--version", action="store", type=int, dest="version", default=1,
            help="The version of the MCM run. Simply an arbitrary integer to allow for keeping track of multiple MCM runs, if desired")
    parser.add_argument("--comment", action="store", dest="comment", default="",
            help="An optional comment to add to the mcmRun DB record")
    args = parser.parse_args()

    mcm_setup(args.letter_code, args.version, args.comment)
