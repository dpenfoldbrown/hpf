
import os
import MySQLdb
import unittest

from hpf.function.createdbs.mygo import Mygo
from hpf.function.term import Term
from hpf.function.createdbs.astral import SuperfamilyEntry
from hpf.function.createdbs.blast_filter import AstralBlastFilter
from hpf.function.createdbs.cdhit_go import CDHitGO
import hpf.function.metric
from hpf.function.predictor import Predictor, get_not_id
import gc


class SimpleFrequencyRunTestCase(unittest.TestCase):
	def tearDown(self):
		print "in tearDown"
		del self.freq_metric
		del self.ba
		del self.records
		del self.filtered
		del self.cd
		del self.clstr_preds
		del self.mg
		gc.collect()
		
	def setUp(self):
		self.freq_metric = hpf.function.metric.Frequency()

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
		self.clstr_preds = self.cd.getClusters(self.filtered, self.mg)

		self.freq_metric.compute_metric(self.clstr_preds)

		frequency_table_name = "test_frequency"
		frequency_single_table_name = "test_frequency_single"

		self.conn = MySQLdb.connect(host=sql_host, user=sql_user, passwd=sql_password, db="hpf")
		self.freq_metric.upload_metric(self.conn, frequency_table_name, delete_table=True)
	

	def testFreqTerm(self):
		assert 4 == self.freq_metric.get_metric(Term(a="all"), Term(a="all")), "wrong frequency count"
		#kdrew: testing get metric by string
		assert 4 == self.freq_metric.get_metric("all", "all"), "wrong frequency count"
		assert 0 == self.freq_metric.get_metric("GO:0000133", "GO:0000135"), "wrong frequency count"
		assert 3 == self.freq_metric.get_metric(Term(a="GO:0004129")), "wrong frequency count"
		assert 3 == self.freq_metric.get_metric("GO:0004129"), "wrong frequency count"
	def testFreqSF(self):
		sf_pred = SuperfamilyEntry(sf_id = "a.51.1")
		assert 3 == self.freq_metric.get_metric(sf_pred), "wrong frequency count"
		assert 3 == self.freq_metric.get_metric(sf_pred, Term(a="all")), "wrong frequency count"
		assert 0 == self.freq_metric.get_metric(sf_pred, Term(a="GO:0001530")), "wrong frequency count"
		assert 2 == self.freq_metric.get_metric(sf_pred, Term(a="GO:0043227")), "wrong frequency count"
		assert 1 == self.freq_metric.get_metric(sf_pred, Predictor(Term(a="GO:0043227").get_not_id())), "wrong frequency count"
		assert 0 == self.freq_metric.get_metric(Predictor(sf_pred.get_not_id()), Term(a="GO:0043227")), "wrong frequency count"
		assert 1 == self.freq_metric.get_metric(Predictor(sf_pred.get_not_id()), Predictor(Term(a="GO:0043227").get_not_id())), "wrong frequency count"
		#kdrew: testing get metric by string
		assert 1 == self.freq_metric.get_metric(get_not_id("a.51.1"), get_not_id("GO:0043227")), "wrong frequency count"
		assert 1 == self.freq_metric.get_metric(get_not_id("a.51.1"), get_not_id("a.51.1")), "wrong frequency count"
		assert 0 == self.freq_metric.get_metric("a.51.1", get_not_id("a.51.1")), "wrong frequency count"


if __name__ == "__main__":
            unittest.main()
    


