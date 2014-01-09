
import getopt
import sys
import logging

def usage():
	print """
		config=(configuration file),
		help(h)=(this message)

		examples:
		python pair_fitter.py --config parameters/test_pair_fitter.py
		"""

def main():
	config_file = None

	try:
		opts, args = getopt.getopt(sys.argv[1:], "hc:", ["help","config="])
	except getopt.GetoptError, exc:
		print "getoptError: ", exc.opt, exc.msg
		usage()
		sys.exit(2)

	for opt, arg in opts:
		if opt in ("--help","-h"):
			usage()
			sys.exit()
		elif opt in ("--config","-c"):
			config_file = arg
		else:
			assert False, logging.error("unhandled option: %s" % opt)

	if config_file != None:
		#kdrew: import config parameters, see florophe_ortho_config.py for format
		logging.debug(config_file)
		#kdrew: module names do not have the .py 
		module_name = config_file.replace(".py","").replace('/','.')
		exec("import "+module_name+" as config")
	else:
		print "config file required"
		usage()
		sys.exit(2)

	
	results_list = []
	for structA in config.structA_list:
		for structB in config.structB_list:

			import hpf.pdb.pairfit as pf
			pfit = pf.PairFit( structA["apdb_filename"], structB["bpdb_filename"], 
								structA["structA_vector_ids"], structB["structB_vector_ids"], 
								selection_id_type = config.selection_id_type, 
								a_atom_names = structA.get('a_atom_names', ['CA','CB'] ), b_atom_names = structB.get('b_atom_names', ['CA','CB'] ), 
								save_dir = config.save_dir )	

			try:
				results = pfit.pair_fit(config.r_length)
			except AttributeError:
				results = pfit.pair_fit()
			#print results
			results_list.append(( structA, structB, results ))

	#kdrew: index [2][0][2] is rmsd, probably should make a dictionary or object of sorts
	results_list.sort( key = lambda k: k[2][0][2] )
	for result in results_list:
		print result[2][0], result[0]["apdb_filename"], result[1]["bpdb_filename"] 

if __name__ == "__main__":
	main()


