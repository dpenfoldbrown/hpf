
import sys
import MySQLdb
import getopt
import ConfigParser
import datetime

from hpf.function.createdbs.mygo import Mygo
from hpf.function.domain import NewAnnotatedDomains, NewAnnotatedDomain
from hpf.function.term import Terms
from hpf.function.prediction import Predictions


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

	fpa = FuncPred_Analysis(config_file, config_section)
	fpa.funcpred_analysis()

class FuncPred_Analysis():

	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		#kdrew: database parameters
		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')
		self.new_mygo_host= config.get(section, 'new_mygo_host')
		self.new_mygo_db = config.get(section, 'new_mygo_db')
		self.new_annotations_table = config.get(section, 'new_annotations_table')
		self.old_mygo_host= config.get(section, 'old_mygo_host')
		self.old_mygo_db = config.get(section, 'old_mygo_db')
		self.old_annotations_table = config.get(section, 'old_annotations_table')

		self.sequence_host = config.get(section, 'sequence_host')
		self.sequence_db = config.get(section, 'sequence_db')
		self.sequence_table = config.get(section, 'sequence_table')
		self.domain_table = config.get(section, 'domain_table')

		self.prediction_host = config.get(section, 'prediction_host')
		self.prediction_db = config.get(section, 'prediction_db')
		self.prediction_table = config.get(section, 'prediction_table')

		#kdrew: mygo parameters
		self.evidence_codes_str = config.get(section,'evidence_codes')
		self.evidence_codes = list(self.evidence_codes_str.split(','))
		self.term_types_str = config.get(section,'term_types')
		self.term_types = list(self.term_types_str.split(','))

		self.pls_min = config.getfloat(section,'pls_min')
		self.base_max = config.getfloat(section,'base_max')
		self.single_domain = config.getboolean(section,'single_domain')

		self.conn = MySQLdb.connect(host=self.sequence_host, user=self.sql_user, passwd=self.sql_password, db=self.sequence_db)
		self.pred_conn = MySQLdb.connect(host=self.prediction_host, user=self.sql_user, passwd=self.sql_password, db=self.prediction_db)
		self.new_mg = Mygo(self.sql_user,self.sql_password, self.new_mygo_host, e_codes = self.evidence_codes, db_name=self.new_mygo_db)
		self.new_mg.connect()
		self.old_mg = Mygo(self.sql_user,self.sql_password, self.old_mygo_host, e_codes = self.evidence_codes, db_name=self.old_mygo_db)
		self.old_mg.connect()

	#kdrew: find difference beween terms using set theory
	#def diff_terms(self, recent_terms, previous_terms):
	#	recent_set = set(recent_terms)
	#	previous_set = set(previous_terms)		
	#	return list(recent_set - previous_set)

	def funcpred_analysis(self):
		cursor = self.conn.cursor(MySQLdb.cursors.DictCursor)

		#kdrew: get domains from diff tables (eg. hddb_IEA_mygoLite_062005_062009)
		domains = NewAnnotatedDomains(self.conn, self.sequence_table, self.domain_table, self.evidence_codes)
		domains.load_domains()

		for domain in domains:

			#kdrew: if single domain flag is set and parent and domain sequence keys do not match, move to the next domain
			#kdrew: generally when parent sequence key and domain sequence key are not the same the protein is a single domain protein
			if self.single_domain and domain.get_parent_sequence_key() != domain.get_domain_sequence_key():
				continue

			print domain.get_parent_sequence_key(), domain.get_domain_sequence_key()

		#kdrew: get_old_annotations
			old_anns = Terms(ecode=self.evidence_codes)
			old_anns.load_terms(conn=self.old_mg.get_conn(), psk=domain.get_parent_sequence_key(), go_type=self.term_types, term_table=self.old_annotations_table)
			#print "old annotations: "
			#old_anns.print_terms()
		#kdrew: get_new_annotations
			new_anns = Terms(ecode=self.evidence_codes)
			new_anns.load_terms(self.new_mg.get_conn(), domain.get_parent_sequence_key(), self.term_types, term_table=self.new_annotations_table)
			#print "new annotations: "
			#new_anns.print_terms()
		#kdrew: get_predictions
			preds = Predictions(conn=self.pred_conn, table=self.prediction_table)
			preds.load_predictions(domain.get_parent_sequence_key(), domain.get_domain_sequence_key())

		#kdrew: predictions_trim = predictions - old_annotations
			predictions_trim = Predictions()
			#kdrew: produces a list of every prediction that is not in the old annotations
			predictions_trim.extend(list((x) for x in preds if x.function not in old_anns))
			#print "predictions_trim: ", predictions_trim.filter(pls_llr=0, base_llr=-2)

		#kdrew: new_annotations_trim = new_annotations - old_annotations
			#kdrew: produces a list of only the new annotations that are truely new, not been seen before
			new_annotations_trim = list(set(new_anns) - set(old_anns))
			#print "new annotations trim: ", new_annotations_trim

		#kdrew: predictions_trim & new_annotations_trim
			new_ann_trim_tmp = list((x) for x in new_annotations_trim if predictions_trim.getByAcc(x.acc))
			#print "new_ann_trim_tmp: ", new_ann_trim_tmp
			prediction_new_annotation_intersection = Predictions()
			prediction_new_annotation_intersection.extend(list((x) for x in predictions_trim if x.function in new_annotations_trim))

			if len(prediction_new_annotation_intersection) > 0:
				print "pred & new: " , prediction_new_annotation_intersection.filter(pls_llr=self.pls_min, base_llr=self.base_max)


if __name__ == "__main__":
	main()



