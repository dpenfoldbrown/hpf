
import os
import sys
import math
import MySQLdb
import getopt
import pickle
import gc
import ConfigParser
from socket import gethostname

from hpf.function.createdbs.blast_filter import AstralBlastFilter
from hpf.function.createdbs.cdhit_go import CDHitGO
from hpf.function.createdbs.mygo import Mygo
import hpf.function.metric


CONFIG_TABLE = "configuration"

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Calculate Metric Tables (python)"

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

	calc = CalcTables(config_file, config_section)
	calc.createTables()

class CalcTables():

	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		#kdrew: blast parameters
		self.my_blast_db = config.get(section, 'blast_db')
		self.my_blast_file = config.get(section, 'blast_file')
		self.my_blast_outfile = config.get(section, 'blast_outfile')
		self.my_blast_exe = config.get(section, 'blast_exe')
		self.e_value_threshold = config.getfloat(section, 'blast_e_value_threshold')
		self.length_threshold = config.getfloat(section, 'blast_length_threshold')
		self.processors = config.getint(section, 'blast_processors')

		#kdrew: cd-hit parameters
		self.my_filtered_fasta = config.get(section, 'cd_hit_filtered_fasta')
		self.my_cd_hit_exe = config.get(section, 'cd_hit_exe')
		self.identity_cutoff = config.getfloat(section, 'cd_hit_identity_cutoff') 
		self.length_cutoff = config.getfloat(section, 'cd_hit_length_cutoff')

		#kdrew: database parameters
		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')
		self.mygo_host= config.get(section, 'mygo_host')
		self.mygo_db = config.get(section, 'mygo_db')
		self.store_host= config.get(section, 'store_host')
		self.store_db = config.get(section, 'store_db')

		#kdrew: store table parameters
		self.frequency_table_name = config.get(section, 'frequency_table_name')
		self.probability_table_name = config.get(section, 'probability_table_name')
		self.log_ratio_table_name = config.get(section, 'log_ratio_table_name')
		self.mi_table_name = config.get(section, 'mi_table_name')

		#kdrew: mygo parameters
		self.evidence_codes_str = config.get(section,'evidence_codes')
		self.evidence_codes = list(self.evidence_codes_str.split(','))

		#kdrew: probability parameters
		self.pseudo_count = config.getint(section,'pseudo_count')

		#kdrew: checkpoint flags for where to start script
		self.blast_checkpoint =  config.getboolean(section,'blast_checkpoint')
		self.filter_checkpoint =  config.getboolean(section,'filter_checkpoint')
		self.cd_hit_checkpoint =  config.getboolean(section,'cd_hit_checkpoint')
		self.go_terms_checkpoint = config.getboolean(section,'go_terms_checkpoint')
		self.metric_checkpoint =  config.getboolean(section,'metric_checkpoint')
		self.store_metric_checkpoint = config.getboolean(section,'store_metric_checkpoint')

		#kdrew: define metric objects
		self.freq_metric = hpf.function.metric.Frequency() 
		self.prob_metric = hpf.function.metric.Probability(pc = self.pseudo_count)
		self.log_ratio_metric = hpf.function.metric.LogRatios()
		self.mi_metric = hpf.function.metric.MutualInformation()

		self.conn = MySQLdb.connect(host=self.store_host, user=self.sql_user, passwd=self.sql_password, db=self.store_db)

		cursor = self.conn.cursor(MySQLdb.cursors.DictCursor)
		create_query = """CREATE TABLE IF NOT EXISTS """+CONFIG_TABLE+""" ( id INT(11) NOT NULL auto_increment, 
			host_name VARCHAR(255), 
			config_file VARCHAR(255), 
			config_section VARCHAR(255), 
			blast_db VARCHAR(255), 
			blast_file VARCHAR(255), 
			blast_outfile VARCHAR(255), 
			blast_exe VARCHAR(255), 
			blast_e_value_threshold DOUBLE, 
			blast_length_threshold DOUBLE, 
			blast_processors INT, 
			cd_hit_filtered_fasta VARCHAR(255), 
			cd_hit_exe VARCHAR(255), 
			cd_hit_identity_cutoff DOUBLE,
			cd_hit_length_cutoff DOUBLE,
			sql_user VARCHAR(255), 
			mygo_host VARCHAR(255), 
			mygo_db VARCHAR(255), 
			store_host VARCHAR(255), 
			store_db VARCHAR(255), 
			frequency_table_name  VARCHAR(255), 
			probability_table_name VARCHAR(255),  
			log_ratio_table_name  VARCHAR(255), 
			mi_table_name VARCHAR(255), 
			evidence_codes text, 
			pseudo_count INT,
			blast_checkpoint VARCHAR(10),  
			filter_checkpoint VARCHAR(10),
			cd_hit_checkpoint VARCHAR(10),
			go_terms_checkpoint  VARCHAR(10),
			metric_checkpoint  VARCHAR(10),
			store_metric_checkpoint  VARCHAR(10),
			timestamp timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
			PRIMARY KEY (id))"""

		print create_query
		cursor.execute(create_query)

		insert_query = """INSERT INTO """+ CONFIG_TABLE + """ 
			( host_name, 
			config_file,
			config_section,
			blast_db, 
			blast_file, 
			blast_outfile, 
			blast_exe, 
			blast_e_value_threshold, 
			blast_length_threshold, 
			blast_processors, 
			cd_hit_filtered_fasta, 
			cd_hit_exe, 
			cd_hit_identity_cutoff, 
			cd_hit_length_cutoff, 
			sql_user, 
			mygo_host,
			mygo_db, 
			store_host, 
			store_db, 
			frequency_table_name, 
			probability_table_name, 
			log_ratio_table_name, 
			mi_table_name, 
			evidence_codes, 
			pseudo_count,
			blast_checkpoint,
			filter_checkpoint,
			cd_hit_checkpoint,
			go_terms_checkpoint, 
			metric_checkpoint,
			store_metric_checkpoint) 
			VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
			"""
			
		host_name = gethostname()
		config_upload = [host_name,
			config_file,	
			section,
			self.my_blast_db, 
			self.my_blast_file, 
			self.my_blast_outfile, 
			self.my_blast_exe, 
			self.e_value_threshold, 
			self.length_threshold, 
			self.processors, 
			self.my_filtered_fasta, 
			self.my_cd_hit_exe, 
			self.identity_cutoff, 
			self.length_cutoff, 
			self.sql_user,
			self.mygo_host,
			self.mygo_db,
			self.store_host,
			self.store_db,
			self.frequency_table_name,
			self.probability_table_name,
			self.log_ratio_table_name,
			self.mi_table_name,
			self.evidence_codes_str,
			self.pseudo_count,
			self.blast_checkpoint,
			self.filter_checkpoint,
			self.cd_hit_checkpoint,
			self.go_terms_checkpoint,
			self.metric_checkpoint,
			self.store_metric_checkpoint ]

		print insert_query
		print config_upload
		cursor.execute(insert_query, config_upload)

	def _default_init(self):
		#kdrew: blast parameters
		self.my_blast_db = "/Users/kdrew/astral/1.75/astral95.1.75"
		#self.my_blast_file = "/Users/kdrew/scripts/function_prediction_python/createdbs/data/go_seq.fasta"
		self.my_blast_file = "/Users/kdrew/data/go/mygoLite_062009/goLite_062009_seq.fasta"
		#self.my_blast_outfile = None
		#self.my_blast_outfile = "/Users/kdrew/data/hpf/go_seq.blast.xml"
		self.my_blast_outfile = "/Users/kdrew/data/go/mygoLite_062009/goLite_062009_seq.blast.xml"
		#self.my_blast_exe = "/usr/bin/blastall"
		self.my_blast_exe = "/Users/patrick/.local/share/blast-2.2.18/bin/blastall"
		self.e_value_threshold = 1e-8
		self.length_threshold = .85
		self.processors = 4

		#kdrew: cd-hit parameters
		self.my_filtered_fasta = self.my_blast_file.rpartition('.')[0]+"_filtered.fasta"
		#self.my_fasta_file = "/Users/kdrew/scripts/function_prediction_python/createdbs/data/go_seq_filtered.fasta"
		self.my_cd_hit_exe = "/Users/kdrew/programs/cd-hit/cd-hit -M 1000"
		self.identity_cutoff = .95
		self.length_cutoff = .5

		#kdrew: database parameters
		self.sql_user="kdrew"
		self.sql_password="kdrew_nyu"
		self.mygo_host="localhost"
		self.mygo_db = "mygoLite_062009"
		#self.store_host="carl.bio.nyu.edu"
		self.store_host="localhost"
		self.store_db = "functionTables"

		#kdrew: store table parameters
		self.frequency_table_name = "frequency_goLite_062009"
		#self.frequency_single_table_name = "frequency_single_goLite_062009"
		self.probability_table_name = "probability_goLite_062009"
		self.log_ratio_table_name = "log_ratio_goLite_062009"
		self.mi_table_name = "mutual_info_goLite_062009"

		#kdrew: mygo parameters
		self.evidence_codes = ['TAS','IDA','IMP','EXP','IPI','IGI','IEP','ISS','ISO','ISA','ISM','IGC','RCA','NAS','IC','ND']

		#kdrew: probability parameters
		self.pseudo_count = 2

		#kdrew: checkpoint flags for where to start script
		self.blast_checkpoint = False
		self.filter_checkpoint = True
		self.cd_hit_checkpoint = True
		self.go_terms_checkpoint = True
		self.metric_checkpoint = True
		self.store_metric_checkpoint = True

		#kdrew: define metric objects
		self.freq_metric = hpf.function.metric.Frequency()
		self.prob_metric = hpf.function.metric.Probability(pc = self.pseudo_count)
		self.log_ratio_metric = hpf.function.metric.LogRatios()
		self.mi_metric = hpf.function.metric.MutualInformations()

		#kdrew: pickle files
		self.filter_pickle_filename = "filter.pkl"
		self.cluster_pickle_filename = "cluster.pkl"
		self.frequency_pickle_filename = "frequency.pkl"
		self.probability_pickle_filename = "probability.pkl"
		self.log_ratio_pickle_filename = "log_ratio.pkl"
		self.mutual_info_pickle_filename = "mutual_info.pkl"

		#self.frequency_pickle_infilename = "pickled_objects/frequency.pkl"
		#self.frequency_pickle_infile= open(self.frequency_pickle_infilename, 'rb') 

		self.filter_pickle_file= open(self.filter_pickle_filename, 'wb')
		self.cluster_pickle_file= open(self.cluster_pickle_filename, 'wb') 
		self.frequency_pickle_file= open(self.frequency_pickle_filename, 'wb') 
		self.probability_pickle_file= open(self.probability_pickle_filename, 'wb')
		self.log_ratio_pickle_file= open(self.log_ratio_pickle_filename, 'wb')
		self.mutual_info_pickle_file= open(self.mutual_info_pickle_filename, 'wb')

		self.conn = MySQLdb.connect(host=self.store_host, user=self.sql_user, passwd=self.sql_password, db=self.store_db)

	def createTables(self):

		self.ba = AstralBlastFilter(self.my_blast_exe, self.my_blast_db, self.my_blast_file, self.e_value_threshold, self.length_threshold, blast_processors=self.processors)

		if self.blast_checkpoint:
			self.records = self.ba.runBlast()
		#kdrew: a little kludgey, if not running blast but are going to filter from blast file, read in the blast out file
		elif self.filter_checkpoint:
			outfile_handle = open(self.my_blast_outfile)
			self.records = self.ba.runBlast(result_handle=outfile_handle)

		if self.filter_checkpoint:
			self.filtered = self.ba.filterBlast(self.records)
			#self.filtered = pickle.load(self.filter_pickle_file)
			#pickle.dump(self.filtered,self.filter_pickle_file)
			
			self.ba.writeFasta(self.filtered, self.my_filtered_fasta)
			outfile_handle.close()

			del self.records
			del self.ba
			gc.collect()

		self.cd = CDHitGO(self.my_cd_hit_exe, self.my_filtered_fasta , id_c=self.identity_cutoff, len_c=self.length_cutoff)
		if self.cd_hit_checkpoint:
			#kdrew: run cluster
			self.cd.runCDHit()

		if self.go_terms_checkpoint:
			#kdrew: get cluster terms
			self.mg = Mygo(self.sql_user,self.sql_password, self.mygo_host, e_codes = self.evidence_codes, db_name=self.mygo_db)
			self.mg.connect()
			#kdrew: if freq_metric is given, the frequency tables are automatically updated and the return from getClusters is an empty dictionary
			self.clstr_preds = self.cd.getClusters(self.filtered, self.mg, self.freq_metric)

			#pickle.dump(self.clstr_preds,self.cluster_pickle_file)
			del self.filtered
			del self.mg
			gc.collect()

		if self.metric_checkpoint:
			self.freq_metric.upload_metric(self.conn, self.frequency_table_name, delete_table=True)

			self.prob_metric.compute_metric(self.freq_metric, pseudo_count=True)
			del self.freq_metric
			gc.collect()
			self.prob_metric.upload_metric(self.conn, self.probability_table_name, delete_table=True)

			self.log_ratio_metric.compute_metric(self.prob_metric)
			self.log_ratio_metric.upload_metric(self.conn, self.log_ratio_table_name, delete_table=True)
			del self.log_ratio_metric
			gc.collect()

			self.mi_metric.compute_metric(self.prob_metric)
			self.mi_metric.upload_metric(self.conn, self.mi_table_name, delete_table=True)

	

if __name__ == "__main__":
	main()



