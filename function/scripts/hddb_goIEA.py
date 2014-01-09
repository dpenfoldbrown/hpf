
import sys
import MySQLdb
import getopt
import ConfigParser
import datetime

from hpf.function.createdbs.blast_filter import GOBlastFilter
from hpf.function.createdbs.mygo import Mygo

EXECUTEMANY_SIZE = 50

CONFIG_TABLE = "configuration"

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to IEA Annotation"

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

	gi = GO_IEA(config_file, config_section)
	gi.goIEA()


class GO_IEA():

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
		self.multi_hits = config.getboolean(section, 'multi_hits')

		self.blast_checkpoint =  config.getboolean(section,'blast_checkpoint')
		self.filter_checkpoint =  config.getboolean(section,'filter_checkpoint')

		#kdrew: database parameters
		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')
		self.mygo_host= config.get(section, 'mygo_host')
		self.mygo_db = config.get(section, 'mygo_db')
		self.store_host= config.get(section, 'store_host')
		self.store_db = config.get(section, 'store_db')

		#kdrew: store table parameters
		self.store_table = config.get(section, 'store_table_name')

		self.store_evidence_code = config.get(section, 'store_evidence_code')
		self.store_source = config.get(section, 'store_source')
		
		#kdrew: mygo parameters
		self.evidence_codes_str = config.get(section,'evidence_codes')
		self.evidence_codes = list(self.evidence_codes_str.split(','))

		self.conn = MySQLdb.connect(host=self.store_host, user=self.sql_user, passwd=self.sql_password, db=self.store_db)
		self.mg = Mygo(self.sql_user,self.sql_password, self.mygo_host, e_codes = self.evidence_codes, db_name=self.mygo_db)
		self.mg.connect()

	def goIEA(self):

		self.ba = GOBlastFilter(self.my_blast_exe, self.my_blast_db, self.my_blast_file, self.e_value_threshold, self.length_threshold, blast_processors=self.processors, multi_hits=self.multi_hits)

		if self.blast_checkpoint:
			self.records = self.ba.runBlast()
		#kdrew: a little kludgey, if not running blast but are going to filter from blast file, read in the blast out file
		elif self.filter_checkpoint:
			outfile_handle = open(self.my_blast_outfile)
			self.records = self.ba.runBlast(result_handle=outfile_handle)

		for blast_record in self.records:
			filtered_record = self.ba.filterBlastOne(blast_record)
			if None != filtered_record:
				if self.multi_hits:
					seqs = list()
					for f_rec in filtered_record:
						seqs.append(f_rec.hit_id)
					terms = self.mg.getGoTerms(seqs, full=True, ancestors=False)
					self.store_terms(filtered_record[0].query_id, terms, seqs)
				else:
					print filtered_record.hit_id
					terms = self.mg.getGoTerms(filtered_record.hit_id, full=True, ancestors=False)
					self.store_terms(filtered_record.query_id, terms, [filtered_record.hit_id])


	def store_terms(self, query_id, terms, hit_seqs):
		cursor = self.conn.cursor(MySQLdb.cursors.DictCursor)
		create_query = """CREATE TABLE IF NOT EXISTS """+self.store_table+""" (
		  id int(11) NOT NULL auto_increment,
		  sequence_key int(11) NOT NULL default '0',
		  domain_sequence_key int(11) NOT NULL,
		  acc varchar(20) NOT NULL default '',
		  name varchar(255) NOT NULL,
		  term_type enum('molecular_function','biological_process','cellular_component') NOT NULL,
		  evidence_code varchar(50) NOT NULL default '',
		  xref_dbname varchar(50) NOT NULL,
		  xref_key varchar(50) NOT NULL,
		  probability double NOT NULL,
		  llr double NOT NULL,
		  level int(11) NOT NULL,
		  source varchar(50) NOT NULL default 'unknown',
		  hit_seqs mediumtext NOT NULL default '',
		  insert_date date default NULL,
		  timestamp timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
		  PRIMARY KEY  (id),
		  UNIQUE KEY sequence_key (sequence_key,domain_sequence_key,acc,evidence_code,source),
		  KEY evidence_code (evidence_code)
		)"""
		cursor.execute(create_query)

		query = "INSERT INTO "+self.store_table+" (sequence_key, acc, name, term_type, evidence_code, source,hit_seqs, insert_date) VALUES (%s, %s,%s,%s,%s,%s,%s,%s)"

		query_list = list()
		for term in terms:
			query_list.append((query_id, term.acc, term.name, term.type, self.store_evidence_code, self.store_source,",".join(map(str,hit_seqs)), datetime.date.today()))

		while query_list:
			small_data, query_list = query_list[:EXECUTEMANY_SIZE], query_list[EXECUTEMANY_SIZE:]
			cursor.executemany(query,small_data)
		
class GO_IEA_Just_Store(GO_IEA):
	def __init__(self, conn, store_table, ecode, source):
		self.conn = conn
		self.store_table = store_table
		self.store_evidence_code = ecode
		self.store_source = source
	def goIEA(self):
		print "not implemented"
	

if __name__ == "__main__":
	main()



