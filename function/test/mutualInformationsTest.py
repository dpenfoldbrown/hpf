import os
import MySQLdb
import unittest
import math

from hpf.function.term import Term
from hpf.function.createdbs.blast_filter import AstralBlastFilter
from hpf.function.createdbs.cdhit_go import CDHitGO
from hpf.function.createdbs.mygo import Mygo
import hpf.function.metric
from hpf.function.metric import TINY_NUM

class SimpleMIRunTestCase(unittest.TestCase):
		
	def setUp(self):
		self.freq_metric = hpf.function.metric.Frequency()
		self.prob_metric = hpf.function.metric.Probability()
		self.mi_metric = hpf.function.metric.MutualInformation()

		#kdrew: run blast filter
		self.my_blast_db = "/Users/kdrew/astral/1.75/astral95.1.75"
		self.my_blast_file = "./createdbs/data/test_seq.fasta"
		self.my_blast_exe = "/Users/patrick/.local/share/blast-2.2.18/bin/blastall"
		self.e_value_threshold = 1e-8
		self.length_threshold = .85
		self.ba = AstralBlastFilter(self.my_blast_exe, self.my_blast_db, self.my_blast_file, self.e_value_threshold, self.length_threshold)
		self.records = self.ba.runBlast()
		self.filtered = self.ba.filterBlast(self.records)

		#kdrew: run cluster
		self.my_fasta_file = "./createdbs/data/test_seq_filtered.fasta"
		self.my_cd_hit_exe = "/Users/kdrew/programs/cd-hit/cd-hit"
		self.identity_cutoff = .95
		self.length_cutoff = .5
		self.cd = CDHitGO(self.my_cd_hit_exe, self.my_fasta_file, id_c=self.identity_cutoff, len_c=self.length_cutoff)
		self.cd.runCDHit()

		sql_user="kdrew"
		sql_password="kdrew_nyu"
		sql_host="mcpeepants.bio.nyu.edu"
		self.mg = Mygo(sql_user,sql_password, sql_host, e_codes =['TAS','IDA','IMP','IEA'], db_name="mygo")
		self.mg.connect()
		self.clstr_terms = self.cd.getClusters(self.filtered, self.mg)

		self.freq_metric.compute_metric(self.clstr_terms)
		self.prob_metric.compute_metric(self.freq_metric)
		self.mi_metric.compute_metric(self.prob_metric)

		mi_table_name = "test_mutual_info"
		self.conn = MySQLdb.connect(host=sql_host, user=sql_user, passwd=sql_password, db="hpf")
		self.mi_metric.upload_metric(self.conn, mi_table_name, delete_table=True)


	def testMI(self):
		#kdrew: compute probs ("a.51.1", "GO:0043227") from scratch to calculate mutual information 
		acc_prob = ((3)/(4+TINY_NUM))
		acc2_prob = ((2)/(4+TINY_NUM))
		not_acc_prob = ((1)/(4+TINY_NUM))
		not_acc2_prob = ((2)/(4+TINY_NUM))

		acc_acc2_prob = ((2)/(3+TINY_NUM)) * acc_prob
		not_acc_acc2_prob =((0)/(1+TINY_NUM)) * not_acc_prob
		acc_not_acc2_prob =((1)/(2+TINY_NUM)) * not_acc2_prob

		#kdrew: mutual information calculation
		x = (acc_acc2_prob * math.log(acc_acc2_prob/(acc_prob*acc2_prob+TINY_NUM)+TINY_NUM)) + (not_acc_acc2_prob * math.log(not_acc_acc2_prob/(not_acc_prob*acc2_prob+TINY_NUM)+TINY_NUM))
		#x += (acc_not_acc2_prob * math.log(acc_not_acc2_prob/(acc_prob*not_acc2_prob+TINY_NUM)+TINY_NUM))

		print "MI(a.51.1,GO:0043227) test: ", x, " compute: ", self.mi_metric.get_metric("a.51.1", "GO:0043227")

		assert  round(x,3) == round(self.mi_metric.get_metric("a.51.1", "GO:0043227"),3), "wrong mutual information"

		#kdrew: compute probs ("all", "all") from scratch to calculate mutual information 
		acc_prob = ((4)/(4+TINY_NUM))
		acc2_prob = ((4)/(4+TINY_NUM))
		not_acc_prob = ((0)/(4+TINY_NUM))
		not_acc2_prob = ((0)/(4+TINY_NUM))

		acc_acc2_prob = ((4)/(4+TINY_NUM)) * ((4)/(4+TINY_NUM))
		not_acc_acc2_prob = ((0)/(4+TINY_NUM)) * acc_prob
		acc_not_acc2_prob = ((0)/(4+TINY_NUM)) * acc_prob

		x = (acc_acc2_prob * math.log(acc_acc2_prob/(acc_prob*acc2_prob+TINY_NUM)+TINY_NUM)) + (not_acc_acc2_prob * math.log(not_acc_acc2_prob/(not_acc_prob*acc2_prob+TINY_NUM)+TINY_NUM))
		#x += (acc_not_acc2_prob * math.log(acc_not_acc2_prob/(acc_prob*not_acc2_prob)))

		print "MI(all,all) test: ", x, " compute: ", self.mi_metric.get_metric("all", "all")

		assert  round(x,3) == round(self.mi_metric.get_metric(Term(a="all"), Term(a="all")),3), "wrong mutual information"

	def testTheoryMITerm(self):
		assert round(0.0,3) == round(self.mi_metric.get_metric("all", "all"),3), "wrong mutual information"
		assert round(0.0,3) == round(self.mi_metric.get_metric("GO:0000133", "GO:0000135"),3), "wrong mutual information"

	def testTheoryMISF(self):
		assert round(0.14384,3) == round(self.mi_metric.get_metric("a.51.1", "GO:0043227"),3), "wrong mutual information"
		assert round(0.04247,3) == round(self.mi_metric.get_metric("GO:0043227", "a.51.1"),3), "wrong mutual information"
		assert round(0.21576,3) == round(self.mi_metric.get_metric("a.51.1","a.51.1"),3), "wrong mutual information"


if __name__ == "__main__":
            unittest.main()
    


