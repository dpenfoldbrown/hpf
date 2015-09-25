
import unittest
import os
import sys
import MySQLdb
import getopt
import ConfigParser
import datetime

from Bio.Blast import NCBIXML

from hpf.function.pdb_webservice_tools import get_pdb_date, get_min_date, PDBBlastFilter

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Experiment PDB Check"

	config_file = "config.cfg"
	config_section = "Default"

	try:
		print "before getopt"
		opts, args = getopt.getopt(sys.argv[1:], "hc:s:", ["help", "config=", "section="]) 
		print "after getopt"
		print opts
		print args
	except getopt.GetoptError:	
		usage()
		sys.exit(2)
	for opt, arg in opts:
		print opt
		print arg
		if opt in ("-h","--help"):
			usage()
			sys.exit()
		elif opt in ("-c","--config"):
			config_file = arg
		elif opt in ("-s","--section"):
			config_section = arg
		else:
			assert False, "unhandled option"

	ssc = Experiment_PDB_Check(config_file, config_section)
	ssc.runDateCheck()

	

class Experiment_PDB_Check():
	
	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		#kdrew: blast parameters
		self.e_value_threshold = config.getfloat(section, 'e_value_threshold')
		self.length_threshold = config.getfloat(section, 'length_threshold')
		self.identity_cutoff = config.getfloat(section, 'identity_cutoff')
		self.date_cutoff = datetime.date(2005,01,01)

		self.experiment_key = config.getint(section, 'experiment_key')
		self.my_blast_outfile = config.get(section, 'my_blast_outfile')

		self.sql_db = config.get(section, 'sql_db')
		self.sql_host = config.get(section, 'sql_host')
		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')
		self.sql_port = config.getint(section, 'sql_port')

		self.sql_conn = MySQLdb.connect(host=self.sql_host, user=self.sql_user, passwd=self.sql_password, db=self.sql_db, port=self.sql_port)


	def runDateCheck(self):
		cursor = self.sql_conn.cursor(MySQLdb.cursors.DictCursor)

		#query = """
#
#				#AND ss.ascession_date >= '2005-01-01'
#
#				SELECT distinct md.sequence_key, s.sequence
#				FROM protein as p, sequence as s, domain as d, mcmData as md
#				WHERE p.sequence_key = d.parent_sequence_key and d.domain_sequence_key = md.sequence_key and d.domain_sequence_key = s.id
#					AND p.experiment_key = %s
#			"""
#
#		print query
#		cursor.execute(query, (self.experiment_key,))
#		rows = cursor.fetchall ()
#		
#		for row in rows:
#			if len(row["sequence"])  > 0:
				#pdb_dates = get_pdb_date_Blast(row["sequence"], self.e_value_threshold, self.length_threshold, self.identity_cutoff)
				
		ba = PDBBlastFilter(eval_cutoff = self.e_value_threshold, length_cutoff = self.length_threshold, identity_cutoff = self.identity_cutoff, multi_hits=True)
		outfile_handle = open(self.my_blast_outfile)
		blast_records = NCBIXML.parse(outfile_handle)
		filtered = ba.filterBlast(blast_records)
		
		for key in filtered.keys():
			pdb_dates = get_pdb_date(filtered[key], parse_hit_id=True)

			#print pdb_dates
			min_date = get_min_date(pdb_dates)
			if min_date != None and self.date_cutoff > min_date[0]:
				print "contaminated: ", key,  min_date
			elif min_date != None:
				print "clean: ", key, min_date
			else:
				print "clean (mindate is None): ",key


if __name__ == "__main__":
	main()


