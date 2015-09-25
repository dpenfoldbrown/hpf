
import sys
import MySQLdb
import getopt
import ConfigParser
import datetime
from socket import gethostname

from hpf.function.createdbs.blast_filter import GOBlastFilter
from hpf.function.createdbs.cdhit_go import CDHitGO
from hpf.function.term import get_name


CONFIG_TABLE = "configuration"

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Function Prediction Report"

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

	pr = PredictionReport(config_file, config_section)
	pr.report()

class PredictionReport():

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
		self.mygo_db= config.get(section, 'mygo_db')

		#kdrew: function prediction table parameters
		self.prediction_host = config.get(section, 'prediction_host')
		self.prediction_db = config.get(section, 'prediction_db')
		self.prediction_table_name = config.get(section, 'prediction_table_name')

		#kdrew: hpf table parameters
		self.hpf_host = config.get(section, 'hpf_host')
		self.hpf_db = config.get(section, 'hpf_db')

		self.sequence_key_file = config.get(section, 'sequence_key_file')
		self.organism_key= config.getint(section, 'organism_key')
		self.outfile= config.get(section, 'outfile')
		self.base_llr_max = config.getint(section, 'base_llr_max')
		self.pls_llr_min = config.getint(section, 'pls_llr_min')

		self.hpf_conn = MySQLdb.connect(host=self.hpf_host, user=self.sql_user, passwd=self.sql_password, db=self.hpf_db)
		self.conn = MySQLdb.connect(host=self.prediction_host, user=self.sql_user, passwd=self.sql_password, db=self.prediction_db)
		self.mg_conn = MySQLdb.connect(host=self.mygo_host, user=self.sql_user,passwd=self.sql_password,db=self.mygo_db)

		self.seqs = []

	def format_header(self):
		format_string = ""
		format_string = format_string+"hostname: "+ gethostname()+"\n"
		format_string = format_string+"mygo_host: "+self.mygo_host+"\n"
		format_string = format_string+"mygo_db: "+self.mygo_db+"\n"

		#kdrew: function prediction table parameters
		format_string = format_string+"prediction_host: "+self.prediction_host +"\n"
		format_string = format_string+"prediction_db: "+self.prediction_db +"\n"
		format_string = format_string+"prediction_table_name: "+self.prediction_table_name+"\n"

		format_string = format_string+"sequence_key_file: "+self.sequence_key_file +"\n"
		format_string = format_string+"outfile: "+self.outfile+"\n"
		format_string = format_string+"base_llr_max: "+str(self.base_llr_max) +"\n"
		format_string = format_string+"pls_llr_min: "+str(self.pls_llr_min) +"\n"

		format_string = format_string+"\n"
		format_string = format_string+"parent_sequence_key\t domain_sequence_key\t go_acc\t go_name\t pls_score\t base_score\n"

		return format_string

	def get_sequences_from_file(self):
		#kdrew: open file
		f = open(self.sequence_key_file)
		lines = f.readlines()
		f.close()

		for line in lines:
			strip_line = line.strip()
			if strip_line != "":
				self.seqs.append(strip_line)

	def get_sequences_from_db(self):
		cursor = self.hpf_conn.cursor(MySQLdb.cursors.DictCursor)
		query = """SELECT DISTINCT sequence_key FROM protein as p WHERE p.experiment_key = %s"""
		cursor.execute(query, (self.organism_key,))
		rows = cursor.fetchall()

		for row in rows:
			self.seqs.append(row["sequence_key"])
		

	def report(self):
		cursor = self.conn.cursor(MySQLdb.cursors.DictCursor)

		if self.sequence_key_file != "":
			self.get_sequences_from_file()
		if self.organism_key != -1:
			self.get_sequences_from_db()

		#kdrew: for each id lookup predictions
		#kdrew: format predictions

		fout = open(self.outfile, "w")
		fout.write(self.format_header())

		for seq in self.seqs:
			#kdrew: get all sequence keys with a specific term_type (ie. molecular_function) and a given evidence code
			query_funcs = """SELECT * FROM """+self.prediction_table_name+""" AS fp WHERE fp.parent_sequence_key = %s AND fp.base_llr <= %s AND fp.pls_llr >= %s ORDER BY fp.pls_llr DESC""" % (seq, self.base_llr_max, self.pls_llr_min)
			print query_funcs
			cursor.execute(query_funcs)
			func_rows = cursor.fetchall()

			for func_row in func_rows:
				mf_name = get_name(func_row["mf_acc"], self.mg_conn)
				out_string = str(func_row["parent_sequence_key"])+"\t"+str(func_row["domain_sequence_key"])+ "\t"+func_row["mf_acc"]+ "\t"+mf_name+ "\t"+str(func_row["pls_llr"])+ "\t"+str(func_row["base_llr"])+"\n"
				fout.write(out_string)


if __name__ == "__main__":
	main()



