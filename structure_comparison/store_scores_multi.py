
import os
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
	parser.add_argument("--results_dir", action="store", dest="results_dir", required=True,
			help="The mammoth results directory")
	parser.add_argument("--structure_type", action="store", dest="structure_type", required=True,
			help="Type of structure (astral, pdb, denovo, other)")
	parser.add_argument("--mammoth_table", action="store", dest="table_destination", default="default",
			help="The StructureMammoth table in hpf database to store mammoth results, default=StructureMammoth, 1186=yeast experiment")


	args = parser.parse_args()

	records = session.query(MammothRun).filter(MammothRun.supergroup_key==args.supergroup_key)\
										.filter(MammothRun.version==args.version)\
										.filter(MammothRun.status=='notuploaded')\
										.all()

	for mammoth_record in records:
		results_file = os.path.join(args.results_dir,"mammoth_sg_{0}_g_{1}_{2}.results".format(mammoth_record.supergroup_key, mammoth_record.group_key1, mammoth_record.group_key2))
		print "loading: " + results_file
		ret = store_scores(results_file, mammoth_record, args.structure_type, args.version, args.table_destination)
		
		mammoth_record.status  = 'complete'
		mammoth_record.comment = "Mammoth run successful, imported results using store_scores.py"
		session.flush()




