#
# A tool to organize the run functionality of Mammoth, corresponding
# to the HPF DB system to keep track of Mammoth runs and avoid redundancy.
# 
# @auth: kdrew 3/24/2013
#
# NOTE: does a locking fetch of what code to work on from the HPF DB (table mammothRun).
# This is cool, as it avoids the wild task pool and PBS/SGE/etc array job and env
# var setup we were previously using. Should set ginzu up this way.
#

import os
import re
import sys
import argparse
import traceback
from hpf.mcm.mammoth import MastodonCL, MammothSimpleParser
from hpf.hddb.db import Session, StructureMammoth, StructureMammoth1186, push_to_db, StructureMammoth_class_lookup
session = Session()


#DBSTORE = True	  # cmdline parameter?
CLEANUP = True	  # cmdline parameter?
DEBUG   = True	  # cmdline parameter?


def main(supergroup_key, version, host, structures_dir, work_base, results_dir, structure_type, retry, dbstore=False, table_destination="default", executable="mammoth"):
	"""
	Organizes the main functionality of getting a MammothRun to work on
	Everything after fetching the mammothRun DB
	record is in a try-catch block, to handle setting DB record status as well as possible.
	"""
	from datetime import datetime
	from hpf.pdb.mammoth import Mammoth

	# Get MammothRun record to work on
	print "DRIVER: Fetching MammothRun record from DB..."
	mammoth_record = fetch_mammothRun_record(session, supergroup_key, version, host, retry)
	if not mammoth_record:
		raise Exception("No MammothRun record supergroup_key {0}, version {1} available to run (batch may be complete!)".format(supergroup_key, version))
	print "DRIVER: Obtained MammothRun id {0} for supergroup_key: {1}, groups: {2},{3}".format(mammoth_record.id, mammoth_record.supergroup_key, mammoth_record.group_key1, mammoth_record.group_key2)
	
	# Set up environment: make working directory and change there
	try:
		work_dir = work_base
		#work_dir = os.path.join(work_base, mammoth_record.prediction_code)
		#if not os.path.isdir(work_dir): os.mkdir(work_dir)
		os.chdir(work_dir)
		print "DRIVER: Working directory holding all temp files: {0}".format(work_dir)
	except Exception as e:
		set_error(session, mammoth_record, e, error_str="Error creating environment")
		raise e

	#kdrew: run mammoth
	output_file=os.path.join(results_dir,"mammoth_sg_{0}_g_{1}_{2}.results".format(mammoth_record.supergroup_key, mammoth_record.group_key1, mammoth_record.group_key2))
	# Create mastodon (mammoth w/contact order) command-line object and run Mammoth with it
	mastodon_cl = MastodonCL(contact_order=True, 
							 experiment=os.path.join(structures_dir,str(mammoth_record.group_key1),"mammoth.list"),
							 prediction=os.path.join(structures_dir,str(mammoth_record.group_key2),"mammoth.list"), 
							 cwd=os.getcwd(), 
							 output=output_file,
							 command=executable, 
							 level=0, verbose=0)
	#DEBUG
	if DEBUG: print "Mammoth command: ", mastodon_cl
	ret = Mammoth(mastodon_cl, parse=False).run()
	
	# On successful finish, set finished time and status (only complete if results in DB) and exit
	if "NORMAL_EXIT" not in open(output_file,'r').read():
		mammoth_record.status  = 'error'
		mammoth_record.comment = "No NORMAL_EXIT, failed"
	elif dbstore:
		store_scores(output_file,mammoth_record, structure_type, version, table_destination)
		mammoth_record.status  = 'complete'
		mammoth_record.comment = "Mammoth run successful"
	else: 
		mammoth_record.status  = 'notuploaded'
		mammoth_record.comment = "Mammoth run successful, results not saved to DB"
	session.flush()

	print "DRIVER: MammothRun id {0} for supergroup_key: {1}, groups: {2},{3} complete".format(mammoth_record.id, mammoth_record.supergroup_key, mammoth_record.group_key1, mammoth_record.group_key2)


def store_scores(output_file, mammoth_record, structure_type, version, table_destination="default"):

    if not os.path.isfile(output_file):
        print "Could not parse Mammoth results: no result file '{0}' (try running first)".format(output_file)
        set_error(session, mammoth_record, e, error_str="Error running Mammoth")
        return None
    with open(output_file) as out_handle:
        scores = list(MammothSimpleParser().parse(out_handle, debug=DEBUG))

    insert_list = []
    for score in scores:                         
        prediction_id = int(score.prediction.split('.')[0])
        experiment_id = int(score.experiment.split('.')[0])

        if prediction_id > experiment_id:
            tmp = prediction_id
            prediction_id = experiment_id
            experiment_id = tmp

        #kdrew: do not check, just try and store
        ##kdrew: check to see if score is already in db
        #sm_pe = session.query(StructureMammoth).filter((StructureMammoth.prediction_id==prediction_id) & (StructureMammoth.experiment_id==experiment_id)).filter(StructureMammoth.version==version).first()
        ## Flipped
        #sm_ep = session.query(StructureMammoth).filter((StructureMammoth.prediction_id==experiment_id) & (StructureMammoth.experiment_id==prediction_id)).filter(StructureMammoth.version==version).first()
        ##session.close()

        #if sm_pe or sm_ep:
        #    print "Score pair already stored in db, skipping {0}, {1}".format(prediction_id, experiment_id)
        #else:
            
        #sm_dbo = StructureMammoth(prediction_id=prediction_id,
        #                          experiment_id=experiment_id,
        #                          prediction_type=structure_type,
        #                          experiment_type=structure_type,
        #                          ini_psi=score.psi1,
        #                          end_psi=score.psi2,
        #                          zscore=score.Zscore,
        #                          evalue=score.Evalue,
        #                          version=version
        #                         )
        #print sm_dbo
        #push_to_db(session, sm_dbo, raise_on_duplicate=False)

        insert_list.append({"prediction_id":prediction_id,
                            "experiment_id":experiment_id,
                            "prediction_type":structure_type,
                            "experiment_type":structure_type,
                            "ini_psi":score.psi1,
                            "end_psi":score.psi2,
                            "zscore":score.Zscore,
                            "evalue":score.Evalue,
                            "version":version})



	#kdrew: lookup class name to insert into desired table
	SM = StructureMammoth_class_lookup[table_destination]
    session.execute(SM.__table__.insert().prefix_with("IGNORE"), insert_list)
	

def set_error(session, record, exception, error_str="Error"):
	"""Sets the error in the MammothRun DB record and prints error messages"""
	record.status  = 'error'
	record.comment = "{0}: {1}".format(error_str, exception)
	session.flush()
	print "DRIVER: MammothRun id {0} for supergroup_key: {1}, groups: {2},{3} failed".format(record.id, record.supergroup_key, record.group_key1, record.group_key2)
	print "DRIVER: Setting MammothRun status to error, halting with exception:\n\n"


def fetch_mammothRun_record(session, supergroup_key, version, host, redo_error):
	"""
	First attempts to lock the hpf.mammothRun table in order to ensure this process gets
	a mammoth record that no other mammothrun_driver process is working on simultaneously.
	
	When a lock is acheived (it's a blocking call), queries the mammothRun table for records
	of the supergroup_key and version specified with a status of 'unprocessed' (if 
	redo_error is True, then also with status 'error').

	If no record, returns None. Otherwise, sets the record's status to 'running',
	sets the host field, unlocks the table, and returns the record.
	"""
	from sqlalchemy import or_
	from hpf.hddb.db import MammothRun
	from sqlalchemy.sql.expression import func
	session.execute("LOCK TABLES {0} WRITE".format(MammothRun.__tablename__))
	if redo_error:
		ready_record = session.query(MammothRun).filter(MammothRun.supergroup_key==supergroup_key)\
											.filter(MammothRun.version==version)\
											.filter(or_(MammothRun.status=='unprocessed', MammothRun.status=='error'))\
											.order_by(func.rand())\
											.first()
	else:
		ready_record = session.query(MammothRun).filter(MammothRun.supergroup_key==supergroup_key).filter(MammothRun.version==version)\
											.filter(MammothRun.status=='unprocessed')\
											.first()
	if ready_record:
		ready_record.status = 'running'
		ready_record.host = host
		session.flush()
	session.execute("UNLOCK TABLES")
	return ready_record


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Driver tool to get MammothRun record from DB (see mammothrun_setup) and run Mammoth on it")
	parser.add_argument("--supergroup_key", action="store", type=int, dest="supergroup_key", required=True,
			help="The supergroup_key to mammoth on. ")
	parser.add_argument("--version", action="store", type=int, dest="version", default=1,
			help="The version of the MCM run. Simply an arbitrary integer to allow for keeping track of multiple MCM runs, if desired")
	parser.add_argument("--host", action="store", dest="host", required=True,
			help="A string to identify where the driver is being run. Machine name, cluster name, etc.")
	parser.add_argument("--work_dir", action="store", dest="work_dir", default="/tmp/mammoth",
			help="The base directory to create a working Mammoth directory (workdir/<code>) in")
	parser.add_argument("--results_dir", action="store", dest="results_dir", required=True,
			help="The directory containing the gzipped rosetta prediction files, eg: ox999.result.gz")
	parser.add_argument("--structures_dir", action="store", dest="structures_dir", required=True,
			help="The directory containing the structures (directory structure can be created by get_structures.py)")
	parser.add_argument("--structure_type", action="store", dest="structure_type", required=True,
			help="Type of structure (astral, pdb, denovo, other)")
	parser.add_argument("--redo_error", action="store_true", dest="redo_error", required=False,
			help="If present, will cause driver to run Mammoth on MammothRun DB entries with 'error' status")
	parser.add_argument("--dbstore", action="store_true", dest="dbstore", required=False,
			help="If present, will upload to database")
	parser.add_argument("--mammoth_table", action="store", dest="table_destination", default="default",
			help="The StructureMammoth table in hpf database to store mammoth results, default=StructureMammoth, 1186=yeast experiment")

	args = parser.parse_args()

	# Set and check directories
	args.work_dir = os.path.abspath(os.path.expanduser(args.work_dir))
	args.results_dir = os.path.abspath(os.path.expanduser(args.results_dir))
	if not os.path.isdir(args.work_dir):
		raise Exception("Working directory {0} is not valid".format(args.work_dir))
	if not os.path.isdir(args.results_dir):
		raise Exception("Results directory {0} is not valid".format(args.results_dir))

	main(args.supergroup_key, args.version, args.host, args.structures_dir, args.work_dir, args.results_dir, args.structure_type, args.redo_error, args.dbstore, args.table_destination)

