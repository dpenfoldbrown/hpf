#
# A setup tool to initialize the HPF DB for an Mammoth run. 
# Adds and inits entries for structure_key in supergroup to the hpf.mammothRun and hpf.mammothGroup tables,
#
# @auth: kdrew 3/23/2013
#
# NOTE: version (here and in the HPF DB) could be overwritten to support a more
# robust and legit version number/system, such as that for the ginzuRun ginzu_version
#

import re
import argparse
import itertools as it


def mammothrun_setup(structure_key_file, supergroup_name, supergroup_comment,  version, comment):
	"""
	"""
	from datetime import datetime
	from hpf.hddb.db import Session, MammothRun, MammothGroup, MammothSupergroup, push_to_db
	session = Session()

	# Check for pre-existing McmRun records with same code and version
	#if session.query(MammothRun).filter(MammothRun.supergroup_key==supergroup_key).filter(MammothRun.version==version).first():
	#	raise Exception("Supergroup {0}, version {1} already exists in DB. Use a different version or clean DB").format(supergroup_key, version)

	msg_dbo = MammothSupergroup(name=supergroup_name, comment=supergroup_comment)
	push_to_db(session, msg_dbo, raise_on_duplicate=True)
	print msg_dbo

	supergroup_key = msg_dbo.id


	struct_keys = []
	handle = open(structure_key_file)
	for line in handle.readlines():
		try:
			int(line.rstrip())
			struct_keys.append(line.rstrip())
		except:
			continue

	#seen = set()
    #seen_add = seen.add
	#return [ x for x in seq if x not in seen and not seen_add(x)]
	group_keys = set()
	for structure_key in struct_keys:
		#kdrew: take last two digits of structure_key and that is its group_key
		group_key = str(structure_key)[-2:]
		group_keys.add(group_key)

		#kdrew: insert into MammothGroup
		mg_dbo = MammothGroup(supergroup_key=supergroup_key,
								group_key=group_key,
								structure_key=structure_key
							 )   
		print mg_dbo
		push_to_db(session, mg_dbo, raise_on_duplicate=False)

	for group_key in group_keys:
		mrs_dbo = MammothRun(supergroup_key=supergroup_key,
								group_key1=group_key,
								group_key2=group_key,
								version=version,
								comment=comment
							)
		print mrs_dbo
		push_to_db(session, mrs_dbo, raise_on_duplicate=False)

	for (group_key1, group_key2) in it.combinations(group_keys,2):
		mr_dbo = MammothRun(supergroup_key=supergroup_key,
								group_key1=group_key1,
								group_key2=group_key2,
								version=version,
								comment=comment
							)
		print mr_dbo
		push_to_db(session, mr_dbo, raise_on_duplicate=False)



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Setup tool to initialize an MCM run in the HPF DB")
	parser.add_argument("--structure_key_file", action="store", dest="structure_key_file", required=True,
			help="File of structure_keys which are to be mammothed all vs all")
	#parser.add_argument("--supergroup_key", action="store", dest="supergroup_key", required=True,
	#		help="The supergroup key which uniquely identifies a set of structures which are to be mammothed all vs all")
	parser.add_argument("--supergroup_name", action="store", dest="supergroup_name", 
			help="The name of the supergroup.")
	parser.add_argument("--supergroup_comment", action="store", dest="supergroup_comment", 
			help="The comment for the supergroup.")
	parser.add_argument("--version", action="store", type=int, dest="version", default=1,
			help="The version of the Mammoth run. Simply an arbitrary integer to allow for keeping track of multiple mammoth runs, if desired")
	parser.add_argument("--comment", action="store", dest="comment", default="",
			help="An optional comment to add to the mammothRun DB record")
	args = parser.parse_args()

	mammothrun_setup(args.structure_key_file, args.supergroup_name, args.supergroup_comment, args.version, args.comment)
