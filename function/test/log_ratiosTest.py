
import os
import MySQLdb
import unittest
import math

from hpf.function.term import Term
from hpf.function.createdbs.astral import SuperfamilyEntry
from hpf.function.createdbs.blast_filter import AstralBlastFilter
from hpf.function.createdbs.cdhit_go import CDHitGO
from hpf.function.createdbs.mygo import Mygo
import hpf.function.metric
from hpf.function.metric import TINY_NUM
from hpf.function.predictor import get_not_id

class SimpleLogRatioRunTestCase(unittest.TestCase):
		
	def setUp(self):
		self.freq_metric = hpf.function.metric.Frequency()
		self.prob_metric = hpf.function.metric.Probability()
		self.log_ratio_metric = hpf.function.metric.LogRatios()

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
		del self.clstr_terms
		self.prob_metric.compute_metric(self.freq_metric)
		self.log_ratio_metric.compute_metric(self.prob_metric)

		log_ratio_table_name = "test_log_ratio"

		self.conn = MySQLdb.connect(host=sql_host, user=sql_user, passwd=sql_password, db="hpf")
		self.log_ratio_metric.upload_metric(self.conn, log_ratio_table_name, delete_table=True)

	

	def testLogRatioTerm(self):
		#kdrew: first TINY_NUM is so don't divide by 0 when computing probability
		#kdrew: second TINY_NUM is so don't divide by 0 when computing log ratio
		#kdrew: third TINY_NUM is so don't take log(0)
		print "lr(all|all): " , self.log_ratio_metric.get_metric(Term(a="all"), Term(a="all"))
		print "lr(all|all) compute: ", math.log(((4)/(4+TINY_NUM))/(TINY_NUM+TINY_NUM)+TINY_NUM)
		assert math.log(((4)/(4+TINY_NUM))/(TINY_NUM+TINY_NUM)+TINY_NUM) == self.log_ratio_metric.get_metric(Term(a="all"), Term(a="all")), "wrong log ratio"
		assert math.log((TINY_NUM)) == self.log_ratio_metric.get_metric(Term(a="GO:0000133"), Term(a="GO:0000135")), "wrong log ratio"

	def testLogRatioSF(self):
		print "log_ratio (GO:0043227 | a.51.1): ", self.log_ratio_metric.get_metric(SuperfamilyEntry(sf_id = "a.51.1"), Term(a="GO:0043227"))
		print "compute log_ratio (GO:0043227 | a.51.1): ", math.log(((2+TINY_NUM)/(3+TINY_NUM))/(TINY_NUM))
		assert math.log((((2)/(3+TINY_NUM))/(TINY_NUM+TINY_NUM))+TINY_NUM) == self.log_ratio_metric.get_metric(SuperfamilyEntry(sf_id = "a.51.1"), Term(a="GO:0043227")), "wrong log ratio "
		assert math.log(TINY_NUM) == self.log_ratio_metric.get_metric(SuperfamilyEntry(sf_id = "a.51.1"), Term(a="GO:0001530")), "wrong log_ratio count"
		assert math.log(TINY_NUM) == self.log_ratio_metric.get_metric(Term(a="GO:0001530"), SuperfamilyEntry(sf_id = "a.51.1")), "wrong log_ratio count"

	def testTheoryLRTerm(self):
		assert round(15.4249,3) == round(self.log_ratio_metric.get_metric("all", "all"),3), "wrong log ratio"
		assert round(14.0386,3) == round(self.log_ratio_metric.get_metric("all", "GO:0043170"),3) , "wrong log ratio"
		assert round(0.0,3) == round(self.log_ratio_metric.get_metric("GO:0043170","all"),3) , "wrong log ratio"
		assert round(0.4054,3) == round(self.log_ratio_metric.get_metric("GO:0043170", "GO:0008324"),3), "wrong log ratio"
		assert round(-16.11809,3) == round(self.log_ratio_metric.get_metric("GO:0000133", "GO:0000135"),3), "wrong log ratio"

	def testTheoryLRSF(self):
		assert round(15.0194,3) == round(self.log_ratio_metric.get_metric("a.51.1", "GO:0043227"),3), "wrong log ratio"
		assert round(0.6931,3) == round(self.log_ratio_metric.get_metric("GO:0043227", "a.51.1"),3), "wrong log ratio"
		assert round(15.4249,3) == round(self.log_ratio_metric.get_metric("a.51.1","a.51.1"),3), "wrong log ratio"

if __name__ == "__main__":
            unittest.main()
    


