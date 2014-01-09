"""
Writes structures in given Mammoth Group to file, using directory-subdirectory file hierarchy
"""

from hpf.hddb.db import Session, Structure, MammothGroup

session = Session()

if __name__ == "__main__":
	import argparse, os, hpf.structure_comparison.structure_mammoth as sm

	parser = argparse.ArgumentParser(description="Get structures from db")

	#parser.add_argument("-f", "--structure_key_file", action="store", dest="structure_key_file", default="",
	#	help="File containing structure keys ")
	parser.add_argument("-s", "--supergroup_key", action="store", type=int, dest="supergroup_key",
		help="Primary key of supergroup")
	parser.add_argument("-d", "--dir", action="store", dest="directory", default=None,
		help="Base directory to put structure files in")

	args = parser.parse_args()

	structure_entries = session.query(MammothGroup).filter(MammothGroup.supergroup_key==args.supergroup_key).all()

	subdir_set = set()

	for structure_entry in structure_entries:

		subdirname = os.path.join(args.directory, "{0}".format(str(structure_entry.group_key)))
		subdir_set.add(subdirname)

		if not os.path.exists(subdirname):
			os.makedirs(subdirname)
		filename = os.path.join(subdirname, "{0}.pdb".format(structure_entry.structure_key))
		print filename
		structure = session.query(Structure).get(structure_entry.structure_key)
		sm.write_struct_file(filename, structure)

	for subdir in subdir_set:
		#kdrew:  make mammoth list files
		subdir_pdbs = os.listdir(subdir)
		mammothlist_filename = os.path.join(subdir,"mammoth.list")
		handle = open(mammothlist_filename, "w")
		handle.write("MAMMOTH List\n")
		handle.write(subdir+"\n")
		for pdb in subdir_pdbs:
			handle.write(pdb+"\n")
		handle.close()






