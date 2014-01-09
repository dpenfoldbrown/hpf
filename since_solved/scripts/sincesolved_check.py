
import unittest
import os
import sys
import MySQLdb
import getopt
import ConfigParser
import datetime

from hpf.function.pdb_webservice_tools import get_pdb_date_Blast, get_min_date

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Since Solved Check"

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

	ssc = SinceSolved_Check(config_file, config_section)
	ssc.runSSCheck()

	

class SinceSolved_Check():
	
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

		self.since_solved_db = config.get(section, 'since_solved_db')
		self.since_solved_host = config.get(section, 'since_solved_host')
		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')

		self.since_solved_conn = MySQLdb.connect(host=self.since_solved_host, user=self.sql_user, passwd=self.sql_password, db=self.since_solved_db)


	def runSSCheck(self):
		cursor = self.since_solved_conn.cursor(MySQLdb.cursors.DictCursor)

		query = """

				SELECT distinct mc.cluster_id AS cluster_id, mc.sequence_key, s.sequence, ss.ascession_date
				FROM since_solved.mcm_cluster mc, since_solved.since_solved ss, hpf.sequence as s
				WHERE mc.cluster_id = ss.cluster_id AND mc.sequence_key = s.id
					AND ss.ascession_date >= '2005-01-01'
			"""

		print query
		cursor.execute(query)
		rows = cursor.fetchall ()
		
		for row in rows:
			if len(row["sequence"])  > 0:
				pdb_dates = get_pdb_date_Blast(row["sequence"], self.e_value_threshold, self.length_threshold, self.identity_cutoff)
				min_date = get_min_date(pdb_dates)
				if min_date != None and self.date_cutoff > min_date[0]:
					print row["sequence_key"], min_date, row["sequence"]


if __name__ == "__main__":
	main()


