
import unittest
import os
from hpf.function.createdbs.mygo import Mygo
from hpf.function.createdbs.cdhit_go import CDHitGO
from hpf.function.createdbs.blast_filter import AstralBlastFilter


class SimpleCDHitRunTestCase(unittest.TestCase):
	def setUp(self):
		print "in cdhit test setup"
		self.my_blast_db = "/Users/kdrew/astral/1.75/astral95.1.75"
		self.my_blast_file = "./createdbs/data/test_seq_cdHit.fasta"
		self.my_blast_outfile = "./createdbs/data/test_seq_cdHit.blast.xml"
		#self.my_blast_exe = "/usr/bin/blastall"
		self.my_blast_exe = "/Users/patrick/.local/share/blast-2.2.18/bin/blastall"
		self.e_value_threshold = 1e-8
		self.length_threshold = .85
		self.processors = 3

		self.my_fasta_file = "./createdbs/data/test_seq_cdHit_filtered.fasta"
		self.my_cd_hit_exe = "/Users/kdrew/programs/cd-hit/cd-hit"
		self.identity_cutoff = .95
		self.length_cutoff = .5

		sql_user="kdrew"
		sql_password="kdrew_nyu"
		sql_host="localhost"


		#kdrew: compare to command line
		#blastall -p "blastp" -d /home/kdrew/astral/1.75/astral95.1.75 -i /home/kdrew/scripts/function_prediction_python/createdbs/data/test_seq.fasta -m 7 -e 1e-8

		self.ba = AstralBlastFilter(self.my_blast_exe, self.my_blast_db, self.my_blast_file, self.e_value_threshold, self.length_threshold, blast_processors=self.processors)
		
		self.records = self.ba.runBlast()
		print "after blast"

		self.filtered = self.ba.filterBlast(self.records)
		print "after filter"

		self.ba.writeFasta(self.filtered)

		self.mg = Mygo(sql_user,sql_password, sql_host, e_codes =['TAS','IDA','IMP'], db_name="mygo")
		self.mg.connect()
		print "conneted to mygo"

		self.mgIEA = Mygo(sql_user,sql_password, sql_host, e_codes =['IEA'], db_name="mygo")
		self.mgIEA.connect()

		#kdrew: compare to command line
		#~/cd-hit/cd-hit/cd-hit -i test_seq_filtered.fasta -o test_seq_filtered.cd_hit.out -c .95 -s .5

		self.cd = CDHitGO(self.my_cd_hit_exe, self.my_fasta_file, id_c=self.identity_cutoff, len_c=self.length_cutoff)
		
		self.cd.runCDHit()

		self.cd.printCDHit()

	def testFastaOut(self):
		assert 1395 == os.path.getsize(self.my_fasta_file.rpartition('.')[0]+".cd_hit.out"), 'wrong file size'

	#kdrew: function deprecated
	#def testGetGoTerms(self):

		#clstr_terms = self.cd.getClusterGoTerms(self.mg)
		#assert 0 == len(clstr_terms['196690']), 'incorrect number of terms returned'
		#assert 21 == len(clstr_terms['1605014']), 'incorrect number of terms returned'

	#	clstr_terms = self.cd.getClusterGoTerms(self.mgIEA)
	#	assert 33 == len(clstr_terms['196690']), 'incorrect number of terms returned'
	#	assert 0 == len(clstr_terms['1605014']), 'incorrect number of terms returned'

	def testGetClusters(self):

		clstr_terms = self.cd.getClusters(self.filtered, self.mg, mustHaveGO=False)
		assert 1 == len(clstr_terms['196690']), 'incorrect number of terms returned'
		assert 22 == len(clstr_terms['1605014']), 'incorrect number of terms returned'

		try:
			x = clstr_terms['3']
			assert False, 'sequence was not removed'
		except KeyError:
			assert True

		clstr_terms = self.cd.getClusters(self.filtered, self.mgIEA, mustHaveGO=False)
		assert 34 == len(clstr_terms['196690']), 'incorrect number of terms returned'
		assert 1 == len(clstr_terms['1605014']), 'incorrect number of terms returned'


if __name__ == "__main__":
            unittest.main()
    


