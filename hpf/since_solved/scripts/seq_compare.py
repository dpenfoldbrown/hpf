
import os,sys
import MySQLdb
import getopt
import ConfigParser

import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.Alphabet import IUPAC

from hpf.alignment_utils import percent_identity

from hpf.hypergeometric import gauss_hypergeom

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Since Solved Sequence Compare"

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

	sc = SequenceCompare(config_file, config_section)

	if sc.cluster_id== -1:
		sc.runSinceSolvedSeqCompare()
	else:
		sc.runSequenceCompareOne(sc.cluster_id)


class SequenceCompare():
	
	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		#kdrew: blast parameters
		self.hpf_db = config.get(section, 'hpf_db')
		self.hpf_host = config.get(section, 'hpf_host')
		self.hpf_port = config.getint(section, 'hpf_port')


		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')

		#self.sequence_key = config.getint(section,'sequence_key')
		self.cluster_id = config.getint(section,'cluster_id')

		self.high_sequence_identity = config.getint(section,'high_sequence_identity')
		self.mcm_value_limit = config.getfloat(section,'mcm_value_limit')
		self.date_cutoff = config.get(section,'date_cutoff')

		self.hpf_conn = MySQLdb.connect(host=self.hpf_host, user=self.sql_user, passwd=self.sql_password, db=self.hpf_db, port=self.hpf_port)

	def runSequenceCompareOne(self, cluster_id):
		cursor = self.hpf_conn.cursor(MySQLdb.cursors.DictCursor)
		query = """
				SELECT mc.cluster_id, mc.sequence_key, s.sequence AS query_sequence, a.astral_id, a.sequence AS astral_sequence, mc.probability  
				FROM since_solved.mcm_cluster AS mc, since_solved.mcm_cluster_max AS mcm, pdb.astral AS a, hpf.sequence AS s 
				WHERE s.id = mc.sequence_key AND mc.experiment_astral_ac = a.astral_id 
				AND mc.probability = mcm.max_prob AND mc.cluster_id = mcm.cluster_id
				AND mc.cluster_id = %s
			"""

		#print query
	
		cursor.execute(query,(cluster_id,))
		row = cursor.fetchone()
		#print row
		if row == None:
			return None
		#kdrew: create biopython sequence objects
		query_seq = SeqRecord(Seq(row['query_sequence'], IUPAC.protein), id=str(row['sequence_key']))
		astral_seq = SeqRecord(Seq(row['astral_sequence'], IUPAC.protein), id=str(row['astral_id']))
		return self.alignment_percent_id(query_seq, astral_seq)

	def alignment_percent_id(self,query_seq, astral_seq):
		#kdrew: align sequences, muscle
		cline = MuscleCommandline(clwstrict=True)
		child = subprocess.Popen(str(cline),stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=(sys.platform!="win32"))
		SeqIO.write((query_seq,astral_seq), child.stdin, "fasta")
		child.stdin.close()
		align = AlignIO.read(child.stdout,"clustal")
		#print align
		#kdrew: percent_identity(alignment)
		align_pi = percent_identity(align)[astral_seq.id]
		return align_pi
		

	def runSinceSolvedSeqCompare(self):
		cursor = self.hpf_conn.cursor(MySQLdb.cursors.DictCursor)

		#query = """
		#		SELECT cluster_id, pred_sf, true_sf
		#		FROM since_solved.true_vs_pred_tophit as tpt
		#	"""
		query = """
				SELECT mc.cluster_id AS cluster_id, SUBSTRING_INDEX(mc.experiment_sccs,'.',3) AS pred_sf, SUBSTRING_INDEX(ss.sccs,'.',3) AS true_sf, mc.probability as top_prob
				FROM since_solved.mcm_cluster mc, since_solved.since_solved ss,
					(SELECT cluster_id, MAX(probability) as max_prob from since_solved.mcm_cluster group by cluster_id) as max_table 
				WHERE mc.cluster_id = ss.cluster_id and max_table.cluster_id = mc.cluster_id and mc.probability = max_table.max_prob 
					AND ss.ascession_date >= DATE(%s)
					AND mc.probability > %s
			"""

		#print query
	
		cursor.execute(query, (self.date_cutoff,self.mcm_value_limit))
		rows = cursor.fetchall ()

		high_id_true = 0
		true_total = 0
		false_total = 0
		sampled = 0

		for row in rows:
			#print row['cluster_id']
			pid = self.runSequenceCompareOne(row['cluster_id'])

			print pid
			if pid == None:
				continue
			if pid >= self.high_sequence_identity:
				sampled +=1
				if row['pred_sf'] == row['true_sf']:
					high_id_true +=1

			print row
			if row['pred_sf'] == row['true_sf']:
				true_total +=1
			else:
				false_total +=1

		print "high_id_true: ", high_id_true, " true_total: ", true_total, " false_total: ", false_total, " sampled: ", sampled
		
		print gauss_hypergeom(high_id_true, true_total, false_total, sampled)


if __name__ == "__main__":
	main()


