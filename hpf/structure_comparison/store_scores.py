
import argparse
from hpf.structure_comparison.mammothrun_driver  import store_scores
from hpf.hddb.db import MammothRun, Session

session = Session()



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Tool to store Mammoth results into DB ")
	parser.add_argument("--version", action="store", type=int, dest="version", default=1,
			help="The version of the MCM run. Simply an arbitrary integer to allow for keeping track of multiple MCM runs, if desired")
	parser.add_argument("--supergroup_key", action="store", type=int, dest="supergroup_key", required=True,
			help="The supergroup of the mammoth run. ")
	parser.add_argument("--group_key1", action="store", type=int, dest="group_key1", required=True,
			help="The group key 1 of the mammoth run. ")
	parser.add_argument("--group_key2", action="store", type=int, dest="group_key2", required=True,
			help="The group key 2 of the mammoth run. ")
	parser.add_argument("--results_file", action="store", dest="results_file", required=True,
			help="The mammoth results file")
	parser.add_argument("--structure_type", action="store", dest="structure_type", required=True,
			help="Type of structure (astral, pdb, denovo, other)")
	parser.add_argument("--mammoth_table", action="store", dest="table_destination", default="default",
			help="The StructureMammoth table in hpf database to store mammoth results, default=StructureMammoth, 1186=yeast experiment")


	args = parser.parse_args()

	record = session.query(MammothRun).filter(MammothRun.supergroup_key==args.supergroup_key)\
										.filter(MammothRun.version==args.version)\
										.filter(MammothRun.group_key1==args.group_key1)\
										.filter(MammothRun.group_key2==args.group_key2)\
										.filter(MammothRun.status=='running')\
										.first()
	store_scores(args.results_file, record, args.structure_type, args.version, args.table_destination)
	record.status  = 'complete'
	record.comment = "Mammoth run successful, imported results using store_scores.py"
	
	session.flush()

