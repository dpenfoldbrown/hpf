

from hpf.hddb.db import Session, MammothGroup, StructureMammoth, StructureMammoth1186, StructureMammoth_class_lookup
from hpf.structure_comparison.structure_mammoth import sm_query
from sqlalchemy import distinct
from sqlalchemy.orm import aliased

session = Session()

if __name__ == "__main__":
	import argparse, os, hpf.structure_comparison.structure_mammoth as sm

	parser = argparse.ArgumentParser(description="Make matrix of mammoth z-scores, another file will be created to store the ordered structure_key listing (.key)")

	parser.add_argument("-s", "--supergroup_key", action="store", type=int, dest="supergroup_key", required=True,
		help="Primary key of supergroup")
	parser.add_argument("-f", "--output_file", action="store", dest="output_file", required=True,
		help="Filename for matrix output")
	parser.add_argument("-t", "--triplet_format", action="store_true", dest="triplet_format", default=False, 
		help="Output in triplet format")
	parser.add_argument("-v", "--version", action="store", type=int, dest="version", default=None,
		help="")
	parser.add_argument("--mammoth_table", action="store", dest="mammoth_table", default="default",
			help="The StructureMammoth table in hpf database to retrieve mammoth results, default=StructureMammoth, 1186=yeast experiment")

	args = parser.parse_args()



	structure_entries = session.query(distinct(MammothGroup.structure_key)).filter(MammothGroup.supergroup_key==args.supergroup_key).all()

	#mg1 = aliased(MammothGroup)
	#mg2 = aliased(MammothGroup)
	
	mat_handle = open(args.output_file,"w")
	key_handle = open(args.output_file+".key","w")


	if args.triplet_format:
		for index, entry in enumerate(structure_entries):

			SM = StructureMammoth_class_lookup[args.mammoth_table]
			#kdrew:  get all structure_key pairs with zscore
			structure_zscores = session.query(SM).filter(SM.prediction_id==entry[0]).filter(SM.version==args.version).all()
			#kdrew:  get reverse
			structure_zscores_rev = session.query(SM).filter(SM.experiment_id==entry[0]).filter(SM.version==args.version).all()

			zscore_dict = dict()
			for structZscore in structure_zscores:
				zscore_dict[structZscore.experiment_id] = structZscore.zscore
			for structZscore in structure_zscores_rev:
				zscore_dict[structZscore.prediction_id] = structZscore.zscore

			key_handle.write(str(entry[0])+"\n")

			for index2, entry2 in enumerate(structure_entries):
				score = zscore_dict[entry2[0]]

				if score:
					mat_handle.write(str(index)+"\t"+str(index2)+"\t"+str(score)+"\n")
				else:
					print "No entry for %s %s" % (entry, entry2)

	#kdrew: regular matrix format
	else:
		for entry in structure_entries:
			key_handle.write(str(entry[0])+"\n")

			for entry2 in structure_entries:
				score = zscore_dict[(entry[0],entry2[0])]

				if not score:
					score = zscore_dict[(entry2[0],entry[0])]

				if score:
					mat_handle.write(str(score.zscore)+"\t")
				else:
					mat_handle.write(str(None)+"\t")

			mat_handle.write("\n")

	mat_handle.close()
	key_handle.close()

