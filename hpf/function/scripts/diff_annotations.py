
import sys
import MySQLdb
import getopt
import ConfigParser
import datetime

from hpf.function.createdbs.blast_filter import GOBlastFilter
from hpf.function.createdbs.cdhit_go import CDHitGO
from hpf.function.createdbs.mygo import Mygo, rows2terms
from hpf.function.scripts.hddb_goIEA import GO_IEA_Just_Store


CONFIG_TABLE = "configuration"

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Function Prediction Analysis"

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

	da = DiffAnnotations(config_file, config_section)
	da.diff()

class DiffAnnotations():

	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		#kdrew: database parameters
		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')
		self.mygo_host= config.get(section, 'mygo_host')
		self.mygo_previous_db = config.get(section, 'mygo_previous_db')
		self.mygo_recent_db = config.get(section, 'mygo_recent_db')

		#kdrew: store table parameters
		self.store_host = config.get(section, 'store_host')
		self.store_db = config.get(section, 'store_db')
		self.store_table_name = config.get(section, 'store_table_name')
		self.store_source = config.get(section, 'store_source')
		self.store_evidence_code = config.get(section, 'store_evidence_code')

		self.annotation_host= config.get(section, 'annotation_host')
		self.annotation_db = config.get(section, 'annotation_db')
		self.annotation_table_recent = config.get(section, 'annotation_table_name_recent')
		self.annotation_table_previous = config.get(section, 'annotation_table_name_previous')

		self.previous_source = config.get(section, 'previous_source')
		self.recent_source = config.get(section, 'recent_source')
		
		#kdrew: mygo parameters
		self.recent_evidence_codes_str = config.get(section,'recent_evidence_codes')
		self.recent_evidence_codes = list(self.recent_evidence_codes_str.split(','))
		self.previous_evidence_codes_str = config.get(section,'previous_evidence_codes')
		self.previous_evidence_codes = list(self.previous_evidence_codes_str.split(','))
		self.term_types_str = config.get(section,'term_types')
		self.term_types = list(self.term_types_str.split(','))

		self.conn = MySQLdb.connect(host=self.annotation_host, user=self.sql_user, passwd=self.sql_password, db=self.annotation_db)
		self.store_conn = MySQLdb.connect(host=self.store_host, user=self.sql_user, passwd=self.sql_password, db=self.store_db)
		self.justStore = GO_IEA_Just_Store(self.store_conn, self.store_table_name, self.store_evidence_code, self.store_source)
		self.mg = Mygo(self.sql_user,self.sql_password, self.mygo_host, e_codes = self.recent_evidence_codes, db_name=self.mygo_recent_db)
		self.mg.connect()

	#kdrew: find difference beween terms using set theory
	def diff_terms(self, recent_terms, previous_terms):
		recent_set = set(recent_terms)
		previous_set = set(previous_terms)		
		return list(recent_set - previous_set)

	def diff(self):
		cursor = self.conn.cursor(MySQLdb.cursors.DictCursor)

		#kdrew: get annotations from recent
		#kdrew: get annotations from previous
		#kdrew: compare

		previous_codes_fmt = "'"+"','".join(self.previous_evidence_codes)+"'"
		print self.previous_evidence_codes 
		print previous_codes_fmt
		recent_codes_fmt = "'"+"','".join(self.recent_evidence_codes)+"'"
		print self.recent_evidence_codes 
		print recent_codes_fmt
		type_fmt = "'"+"','".join(self.term_types)+"'"
		print self.term_types
		print type_fmt

		#kdrew: get all sequence keys with a specific term_type (ie. molecular_function) and a given evidence code
		query_seqs = """select distinct hig.sequence_key from """+self.annotation_table_recent+""" as hig where hig.term_type in (%s) and hig.evidence_code in (%s)""" % (type_fmt, recent_codes_fmt)
		print query_seqs
		cursor.execute(query_seqs)
		seqs = cursor.fetchall()
		print len(seqs)

		#kdrew: for every sequence key
		for seq in seqs:
			sequence_key = seq["sequence_key"]
			#kdrew: get all annotations for the sequence key from the earlier database
			query = """select distinct t2.acc, t2.name, t2.term_type
				from """+self.annotation_table_previous+""" as hig, """+self.mygo_previous_db+""".term as t, """+self.mygo_previous_db+""".term as t2, """+self.mygo_previous_db+""".graph_path as gp 
				where gp.term2_id = t.id and t.acc = hig.acc and t2.id = gp.term1_id and hig.term_type in (%s) and t2.term_type in (%s) and hig.evidence_code in (%s) and hig.sequence_key = %s""" % (type_fmt,type_fmt,previous_codes_fmt,str(sequence_key))
			cursor.execute(query)

			previous_term_rows = cursor.fetchall()
			previous_terms = rows2terms(previous_term_rows, self.previous_evidence_codes,full=True)

			#kdrew: get all annotations for the sequence key from the recent database
			query2 = """select distinct t2.acc, t2.name, t2.term_type
				from """+self.annotation_table_recent+""" as hig, """+self.mygo_recent_db+""".term as t, """+self.mygo_recent_db+""".term as t2, """+self.mygo_recent_db+""".graph_path as gp 
				where gp.term2_id = t.id and t.acc = hig.acc and t2.id = gp.term1_id and hig.term_type in (%s) and t2.term_type in (%s) and hig.evidence_code in (%s) and hig.sequence_key = %s""" % (type_fmt,type_fmt,recent_codes_fmt,str(sequence_key))
			cursor.execute(query2)
			recent_term_rows = cursor.fetchall()
			recent_terms = rows2terms(recent_term_rows, self.recent_evidence_codes,full=True)

			new_terms = self.diff_terms(recent_terms, previous_terms)
			print "sequence_key: ", sequence_key, " new terms: ", new_terms
			self.justStore.store_terms(sequence_key, new_terms)


if __name__ == "__main__":
	main()



