##
## Quick script to set up directory structure and create fasta (and ss2, if available)
## files for rosetta ab initio structure prediction
##

import os
import argparse
from hpf.pdb.psipred import HPFPsipredWrap
from hpf.hddb.db import Session, Domain, FilesystemOutfile as Foldable


BASE_DIR ="/home/dpb3/rosetta"
WORK_DIR ="{0}/work".format(BASE_DIR)
NR_DB    ="/data/nr/nr"


def main(prediction_code):
    print "Running Rosetta ab init setup for code {0}".format(prediction_code)

    # Set up DB and query for foldable and domain records
    print "Opening session and querying DB for domain and foldable"
    session = Session()
    domain = session.query(Domain).filter_by(ibm_prediction_code=prediction_code).first()
    foldable = session.query(Foldable).filter_by(prediction_code=prediction_code).first()
    if not (domain and foldable):
        raise Exception("No domain or foldable record (or both) could be found for code {0}".format(prediction_code))
    
    # Check/make work dir
    if not os.path.isdir(BASE_DIR):
        os.mkdir(WORK_DIR)

    # Check NR
    if not os.path.isfile(NR_DB):   
        raise Exception("Given NR database {0} not valid".format(NR_DB))

    # Make and change to working dir
    working_dir = os.path.join(WORK_DIR, prediction_code)
    os.mkdir(working_dir)
    print "Created working directory {0}".format(working_dir)
    os.chdir(working_dir)

    # Create foldable FASTA file
    fasta_file = os.path.join(working_dir, "{0}.fasta".format(prediction_code))
    with open(fasta_file, 'w') as handle:
        handle.write(">domain:{0}|code:{1}|foldable_seq_key:{2}\n".format(domain.id, foldable.prediction_code, foldable.sequence_key))
        handle.write("{0}\n".format(foldable.sequence.sequence))
    print "Created foldable record fasta file {0}".format(fasta_file)

    # Currently DB does not hold psipred ss2-format files (damn). Run psipred, ignore return pred (only need file)
    print "Running Psipred on foldable sequence via HPFPsipredWrap..."
    HPFPsipredWrap(sequence_key=foldable.sequence_key, nr_db=NR_DB, ginzu_version=4, dir=working_dir, autorun=True, dbstore=True, debug=True)

    # Complete and exit
    print "Setup for Rosetta ab initio complete. Clean up files as necessary"
    print "(NOTE: Psipred ss2 outfile (vformat) is the <seqkey>.psipred.2 file)"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Setup for Rosetta ab init struct pred')
    parser.add_argument("--code", action="store", dest='prediction_code', required=True,
            help='Prediction code (eg: aa123) to set up')
    args = parser.parse_args()
    main(args.prediction_code)
