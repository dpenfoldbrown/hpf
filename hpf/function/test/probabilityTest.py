
import os
import MySQLdb
import unittest

from hpf.function.createdbs.mygo import Mygo
from hpf.function.term import Term
from hpf.function.createdbs.astral import SuperfamilyEntry
from hpf.function.createdbs.blast_filter import AstralBlastFilter
from hpf.function.createdbs.cdhit_go import CDHitGO
import hpf.function.metric
from hpf.function.metric import TINY_NUM
from hpf.function.predictor import get_not_id


class SimpleProbabilityRunTestCase(unittest.TestCase):

	def setUp(self):
		self.freq_metric = hpf.function.metric.Frequency()
		self.prob_metric = hpf.function.metric.Probability()

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

		probability_table_name = "test_probability"

		self.conn = MySQLdb.connect(host=sql_host, user=sql_user, passwd=sql_password, db="hpf")
		self.prob_metric.upload_metric(self.conn, probability_table_name, delete_table=True)

		#print "\n\nfreq_metric:"
		#self.freq_metric.printTables()
		#print "\n\nprob_metric:"
		#self.prob_metric.printTables()
		#print "\n\nbgprob_metric:"
		#print "\n\n"
	

	def testProbTerm(self):
		print "testing probability of terms"
		print "P(all|all) ", self.prob_metric.get_metric("all", "all")
		#kdrew: P( all | all)
		#assert (4+TINY_NUM)/(4+TINY_NUM) == self.prob_metric.get_metric("all", "all"), "wrong probability"
		assert (4)/(4+TINY_NUM) == self.prob_metric.get_metric("all", "all"), "wrong probability"

		#kdrew: P( GO:0043170 | all )
		assert (1)/(4+TINY_NUM) == self.prob_metric.get_metric("all", "GO:0043170" ), "wrong probability"

		#kdrew: P( GO:0043170 | GO:0008324)
		assert (1)/(3+TINY_NUM) == self.prob_metric.get_metric(Term(a="GO:0008324"), Term(a="GO:0043170") ), "wrong probability"

		#kdrew: P( GO:0000135 | GO:0000133 ) testing something not there
		assert TINY_NUM == self.prob_metric.get_metric(Term(a="GO:0000133"), Term(a="GO:0000135")), "wrong probability"

	def testProbSF(self):
		#kdrew: P( GO:0043227 | a.51.1 )
		#assert (2+TINY_NUM)/(3+TINY_NUM) == self.prob_metric.get_metric(SuperfamilyEntry(sf_id = "a.51.1"), Term(a="GO:0043227")), "wrong probability"
		assert (2)/(3+TINY_NUM) == self.prob_metric.get_metric(SuperfamilyEntry(sf_id = "a.51.1"), Term(a="GO:0043227")), "wrong probability"

		#kdrew: P( GO:0001530 | a.51.1 )
		assert TINY_NUM == self.prob_metric.get_metric(SuperfamilyEntry(sf_id = "a.51.1"), Term(a="GO:0001530")), "wrong probability"

		#kdrew: P( all | not a.51.1 )
		print "P(all|not a.51.1): ",self.prob_metric.get_metric(SuperfamilyEntry(SuperfamilyEntry(sf_id = "a.51.1").get_not_id()), Term(a="all"))
		assert round(1.0,3) == round(self.prob_metric.get_metric(SuperfamilyEntry(SuperfamilyEntry(sf_id = "a.51.1").get_not_id()), Term(a="all")),3), "wrong probability"
		#assert (1)/(1+TINY_NUM) == self.prob_metric.get_metric(SuperfamilyEntry(SuperfamilyEntry(sf_id = "a.51.1").get_not_id()), Term(a="all")), "wrong probability"


		#kdrew: P( not a.51.1 | all )
		print "P(nota.51.1|all): ", self.prob_metric.get_metric(Term(a="all"), SuperfamilyEntry(SuperfamilyEntry(sf_id = "a.51.1").get_not_id()))
		assert (1.0 - (3.0/(4+TINY_NUM))) == self.prob_metric.get_metric(Term(a="all"), SuperfamilyEntry(SuperfamilyEntry(sf_id = "a.51.1").get_not_id())), "wrong probability"

		#kdrew: P( a.51.1 | a.51.1 )
		assert (3)/(3+TINY_NUM) == self.prob_metric.get_metric("a.51.1","a.51.1"), "wrong probability"

	def testBGProbSF(self):
		#kdrew: P( a.51.1 )
		print "P(a.51.1): ", self.prob_metric.get_metric("a.51.1",None)
		assert (3)/(4+TINY_NUM) == self.prob_metric.get_metric("a.51.1",None), "wrong probability"
		#kdrew: P( not a.51.1 )
		print "not a.51.1: ", self.prob_metric.get_metric(get_not_id("a.51.1"),None)
		assert 1.0 - (3)/(4+TINY_NUM) == self.prob_metric.get_metric(get_not_id("a.51.1"),None), "wrong probability"

		#kdrew: P( a.51.1 )
		assert (3)/(4+TINY_NUM) == self.prob_metric.get_metric(SuperfamilyEntry(sf_id = "a.51.1"),None), "wrong probability"
		#kdrew: P( not a.51.1 )
		assert 1.0 - (3)/(4+TINY_NUM) == self.prob_metric.get_metric(SuperfamilyEntry(SuperfamilyEntry(sf_id = "a.51.1").get_not_id()),None), "wrong probability"

	def testBGProbMF(self):
		#kdrew: P( GO:0043227 )
		assert (2)/(4+TINY_NUM) == self.prob_metric.get_metric( Term(a="GO:0043227"),None), "wrong probability"
		#kdrew: P( not GO:0043227 )
		assert 1.0 - (2)/(4+TINY_NUM) == self.prob_metric.get_metric( Term(Term(a="GO:0043227").get_not_id()),None), "wrong probability"
		#kdrew: P( GO:0001530 )
		assert TINY_NUM == self.prob_metric.get_metric( Term(a="GO:0001530"),None), "wrong probability"

	def testPseudoCountProb(self):
		pseudo_count_test = 2
		self.prob_metric2 = hpf.function.metric.Probability(pc = pseudo_count_test)
		self.prob_metric2.compute_metric(self.freq_metric, pseudo_count=True)

		x = (((2)/(3+TINY_NUM)) * 2) + (((2)/(4+TINY_NUM)) * pseudo_count_test)
		x = x/(2+pseudo_count_test+TINY_NUM)
		assert x == self.prob_metric2.get_metric(SuperfamilyEntry(sf_id = "a.51.1"), Term(a="GO:0043227")), "wrong probability"
		assert x == self.prob_metric2.get_metric("a.51.1", "GO:0043227"), "wrong probability"

		x = (((4)/(4+TINY_NUM)) * 4) + (((4)/(4+TINY_NUM)) * pseudo_count_test)
		x = x/(4+pseudo_count_test+TINY_NUM)
		assert x == self.prob_metric2.get_metric("all","all"), "wrong probability"

		#kdrew: P( GO:0001530 | a.51.1 )
		assert TINY_NUM == self.prob_metric2.get_metric(SuperfamilyEntry(sf_id = "a.51.1"), Term(a="GO:0001530")), "wrong probability"

	def testPseudoCountProb(self):
		pseudo_count_test = 2
		self.prob_metric2 = hpf.function.metric.Probability(pc = pseudo_count_test)
		self.prob_metric2.compute_metric(self.freq_metric, pseudo_count=True)

		keys = self.prob_metric2.get_all_ids()
		for key in keys:
			assert 0.99 < self.prob_metric2.get_metric(key[0],"all")	
		

	def testTheoryProbTerm(self):
		assert round(1.0,3) == round(self.prob_metric.get_metric("all", "all"),3), "wrong probability"
		assert round(0.0,3) == round(self.prob_metric.get_metric("all", get_not_id("all")),3), "wrong probability"
		print "all: ", self.prob_metric.get_metric("all",None)
		print "all|all: ", self.prob_metric.get_metric("all","all")
		print "all|notall: ", self.prob_metric.get_metric(get_not_id("all"),"all")
		assert round(0.0,3) == round(self.prob_metric.get_metric(get_not_id("all"),"all"),3), "wrong probability"
		assert round(0.25,3) == round(self.prob_metric.get_metric("all", "GO:0043170"),3) , "wrong probability"
		assert round(0.3333,3) == round(self.prob_metric.get_metric("GO:0008324", "GO:0043170"),3), "wrong probability"
		assert round(0.0,3) == round(self.prob_metric.get_metric("GO:0000133", "GO:0000135"),3), "wrong probability"

	def testTheoryProbSF(self):
		assert round(0.66666,3) == round(self.prob_metric.get_metric("a.51.1", "GO:0043227"),3), "wrong probability"
		assert round(0.5,3) == round(self.prob_metric.get_metric(get_not_id("GO:0043227"),"a.51.1"),3), "wrong probability"
		assert round(0.0,3) == round(self.prob_metric.get_metric("a.51.1", "GO:0001530"),3), "wrong probability"
		assert round(1.0,3) == round(self.prob_metric.get_metric(get_not_id("a.51.1"), "all"),3), "wrong probability"
		assert round(0.25,3) == round(self.prob_metric.get_metric("all", get_not_id("a.51.1")),3), "wrong probability"
		assert round(1.0,3) == round(self.prob_metric.get_metric("a.51.1","a.51.1"),3), "wrong probability"

	def testTheoryPseudoCountProb(self):
		pseudo_count_test = 2
		self.prob_metric2 = hpf.function.metric.Probability(pc = pseudo_count_test)
		self.prob_metric2.compute_metric(self.freq_metric, pseudo_count=True)
		#self.freq_metric.printTables()
		probability_table_name = "test_pseudo_probability"
		self.prob_metric2.upload_metric(self.conn, probability_table_name, delete_table=True)

		print "psuedo_count P(GO:0043227|a.51.1): ", self.prob_metric2.get_metric("a.51.1", "GO:0043227")
		assert round(0.5833,3) == round(self.prob_metric2.get_metric("a.51.1", "GO:0043227"),3), "wrong probability"

		assert round(1.0,3) == round(self.prob_metric2.get_metric("all", "all"),3), "wrong probability"

		#kdrew: P( GO:0001530 | a.51.1 )
		assert round(0.0,3) == round(self.prob_metric2.get_metric("a.51.1", "GO:0001530"),3), "wrong probability"


	def testSanity(self):
		pseudo_count_test = 2
		self.prob_metric2 = hpf.function.metric.Probability(pc = pseudo_count_test)
		self.prob_metric2.compute_metric(self.freq_metric, pseudo_count=True)


		for key in self.prob_metric.get_all_ids():
			print key,": ",self.prob_metric2.get_metric(key[0],key[1])
			if 1 <= self.prob_metric2.get_metric(key[0],key[1]):
				print "over one ", key,": ",self.prob_metric2.get_metric(key[0],key[1])
			if 0 >= self.prob_metric2.get_metric(key[0],key[1]):
				print "under zero ", key,": ",self.prob_metric2.get_metric(key[0],key[1])
				print "key0: ", self.prob_metric2.get_metric(key[0])
				print "key1: ", self.prob_metric2.get_metric(key[1])
				print "P(1|not0): ", self.prob_metric.get_metric(key[1])
			assert 1 >= self.prob_metric2.get_metric(key[0],key[1]), "metric larger than 1: "
			assert 0 <= self.prob_metric2.get_metric(key[0],key[1]), "metric smaller than 0: "

if __name__ == "__main__":
            unittest.main()
    


