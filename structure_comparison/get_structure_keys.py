"""
Fetches all the relevant MCM (ie, de Novo for unknown domains) structure keys for a given
experiment.
Outputs struct-key per line to given filename
"""

from hpf.hddb.db import Session, FilesystemOutfile, Protein, McmData
from sqlalchemy import distinct

session = Session()

if __name__ == "__main__":
	import argparse, os, hpf.structure_comparison.structure_mammoth as sm

	parser = argparse.ArgumentParser(description="Get structure keys from db")

	#parser.add_argument("-f", "--structure_key_file", action="store", dest="structure_key_file", default="",
	#	help="File containing structure keys ")
	parser.add_argument("-e", "--experiment_key", action="store", type=int, dest="experiment_key", required=True,
		help="Primary key of experiment")
	parser.add_argument("-f", "--structure_key_file", action="store", dest="structure_key_file", required=True,
		help="File to put structure keys files in")
	parser.add_argument("-c", "--ss_class", action="store", dest="ss_class", default=None,
		help="Secondary structure class, (1=alpha, 2=beta, 3=alpha_beta)")

	args = parser.parse_args()

	if not args.ss_class:
		structure_entries = session.query(distinct(McmData.structure_key)).filter(McmData.sequence_key==FilesystemOutfile.sequence_key).filter(FilesystemOutfile.parent_sequence_key==Protein.sequence_key).filter(Protein.experiment_key==args.experiment_key).all()
	else:
		structure_entries = session.query(distinct(McmData.structure_key)).filter(McmData.sequence_key==FilesystemOutfile.sequence_key).filter(FilesystemOutfile.parent_sequence_key==Protein.sequence_key).filter(Protein.experiment_key==args.experiment_key).filter(McmData.__dict__['class']==args.ss_class).all()

	
	handle = open(args.structure_key_file,"w")

	for entry in structure_entries:
		handle.write(str(entry[0])+"\n")

	handle.close()

