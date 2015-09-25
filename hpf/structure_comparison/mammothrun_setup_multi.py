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

def mammothrun_setup_multi_supergroup(supergroup_key1, supergroup_key2, new_supergroup_name, new_supergroup_comment,  version, comment):
	"""
	"""
	from datetime import datetime
	from hpf.hddb.db import Session, MammothRun, MammothGroup, MammothSupergroup, push_to_db
	session = Session()

	msg_dbo = MammothSupergroup(name=new_supergroup_name, comment=new_supergroup_comment)
	push_to_db(session, msg_dbo, raise_on_duplicate=True)
	print msg_dbo

	new_supergroup_key = msg_dbo.id

	sg1_struct_keys = []
	groups1 = session.query(MammothGroup).filter(MammothGroup.supergroup_key == supergroup_key1).all()
	for entry in  groups1:
		sg1_struct_keys.append(entry.structure_key)

	sg2_struct_keys = []
	groups2 = session.query(MammothGroup).filter(MammothGroup.supergroup_key == supergroup_key2).all()
	for entry in  groups2:
		sg2_struct_keys.append(entry.structure_key)

	sg1_group_keys = set()
	for structure_key in sg1_struct_keys:
		#kdrew: take last two digits of structure_key and that is its group_key
		group_key = str(supergroup_key1)+str(structure_key)[-2:]
		sg1_group_keys.add(group_key)

		#kdrew: insert into MammothGroup
		mg_dbo = MammothGroup(supergroup_key=new_supergroup_key,
								group_key=group_key,
								structure_key=structure_key
							 )   
		print mg_dbo
		push_to_db(session, mg_dbo, raise_on_duplicate=False)

	sg2_group_keys = set()
	for structure_key in sg2_struct_keys:
		#kdrew: take last two digits of structure_key and that is its group_key
		group_key = str(supergroup_key2)+str(structure_key)[-2:]
		sg2_group_keys.add(group_key)

		#kdrew: insert into MammothGroup
		mg_dbo = MammothGroup(supergroup_key=new_supergroup_key,
								group_key=group_key,
								structure_key=structure_key
							 )   
		print mg_dbo
		push_to_db(session, mg_dbo, raise_on_duplicate=False)

	for group_key1 in sg1_group_keys:
		for group_key2 in sg2_group_keys:
			mr_dbo = MammothRun(supergroup_key=new_supergroup_key,
									group_key1=group_key1,
									group_key2=group_key2,
									version=version,
									comment=comment
								)
			print mr_dbo
			push_to_db(session, mr_dbo, raise_on_duplicate=False)




if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Setup tool to initialize an MCM run in the HPF DB")
	parser.add_argument("--supergroup_key1", action="store", dest="supergroup_key1", required=True,
			help="The supergroup key which uniquely identifies a set of structures which are to be mammothed all vs all")
	parser.add_argument("--supergroup_key2", action="store", dest="supergroup_key2", required=True,
			help="The supergroup key which uniquely identifies a set of structures which are to be mammothed all vs all")
	parser.add_argument("--supergroup_name", action="store", dest="supergroup_name", 
			help="The name of the supergroup.")
	parser.add_argument("--supergroup_comment", action="store", dest="supergroup_comment", 
			help="The comment for the supergroup.")
	parser.add_argument("--version", action="store", type=int, dest="version", default=1,
			help="The version of the Mammoth run. Simply an arbitrary integer to allow for keeping track of multiple mammoth runs, if desired")
	parser.add_argument("--comment", action="store", dest="comment", default="",
			help="An optional comment to add to the mammothRun DB record")
	args = parser.parse_args()

	mammothrun_setup_multi_supergroup(args.supergroup_key1, args.supergroup_key2, args.supergroup_name, args.supergroup_comment, args.version, args.comment)

