
import unittest
import os
from hpf.function.createdbs.blast_filter import AstralBlastFilter, GOBlastFilter


class SimpleBlastRunTestCase(unittest.TestCase):
	def setUp(self):
		self.my_blast_db = "/Users/kdrew/astral/1.75/astral95.1.75"
		self.my_blast_file = "./createdbs/data/test_seq.fasta"
		self.my_blast_outfile = "./createdbs/data/test_seq.blast.xml"
		#self.my_blast_exe = "/usr/bin/blastall"
		self.my_blast_exe = "/Users/patrick/.local/share/blast-2.2.18/bin/blastall"
		self.e_value_threshold = 1e-8
		self.length_threshold = .85
		self.processors = 3

		#kdrew: compare to command line
		#blastall -p "blastp" -d /home/kdrew/astral/1.75/astral95.1.75 -i /home/kdrew/scripts/function_prediction_python/createdbs/data/test_seq.fasta -m 7 -e 1e-8

		self.ba = AstralBlastFilter(self.my_blast_exe, self.my_blast_db, self.my_blast_file, self.e_value_threshold, self.length_threshold, blast_processors=self.processors)
		
		self.records = self.ba.runBlast()

		self.filtered = self.ba.filterBlast(self.records)
		
		self.ba.writeFasta(self.filtered)


	def testNumFiltered(self):
		assert 5 == len(self.filtered), 'too many records returned'

	def testKeyFiltered(self):
		assert '3' in self.filtered, 'wrong record filtered'
		assert '1' not in self.filtered, 'wrong record filtered'

	def testSuperfamily(self):
		assert 'a.51.1' == self.filtered[str(1605014)].superfamily, 'wrong superfamily'

	def testFastaOut(self):
		assert 1515 == os.path.getsize(self.my_blast_file.rpartition('.')[0]+"_filtered.fasta"), 'wrong file size'

class ReadInBlastRunTestCase(unittest.TestCase):
	def setUp(self):
		self.my_blast_db = "/Users/kdrew/astral/1.75/astral95.1.75"
		self.my_blast_file = "./createdbs/data/test_seq.fasta"
		self.my_blast_outfile = "./createdbs/data/test_seq.blast.xml"
		#self.my_blast_exe = "/usr/bin/blastall"
		self.my_blast_exe = "/Users/patrick/.local/share/blast-2.2.18/bin/blastall"
		self.e_value_threshold = 1e-8
		self.length_threshold = .85
		self.processors = 3

		#kdrew: compare to command line
		#blastall -p "blastp" -d /home/kdrew/astral/1.75/astral95.1.75 -i /home/kdrew/scripts/function_prediction_python/createdbs/data/test_seq.fasta -m 7 -e 1e-8

		self.ba = AstralBlastFilter(self.my_blast_exe, self.my_blast_db, self.my_blast_file, self.e_value_threshold, self.length_threshold, blast_processors=self.processors)

		self.outfile_handle = open(self.my_blast_outfile)
		
		self.records = self.ba.runBlast(self.outfile_handle)

		self.filtered = self.ba.filterBlast(self.records)
		
		self.ba.writeFasta(self.filtered)


	def testFileOpen(self):
		assert self.outfile_handle != None, 'xml file not open'

	def testNumFiltered(self):
		assert 5 == len(self.filtered), 'too many records returned'

	def testKeyFiltered(self):
		assert '3' in self.filtered, 'wrong record filtered'
		assert '1' not in self.filtered, 'wrong record filtered'

	def testSuperfamily(self):
		assert 'a.51.1' == self.filtered[str(1605014)].superfamily, 'wrong superfamily'

	def testFastaOut(self):
		assert 1515 == os.path.getsize(self.my_blast_file.rpartition('.')[0]+"_filtered.fasta"), 'wrong file size'

class GOBlastRunTestCase(unittest.TestCase):
	def setUp(self):
		self.my_blast_db = None
		self.my_blast_file = None
		self.my_blast_outfile = "./test/data/hddb_test.blast.xml"
		self.my_blast_exe = None
		self.e_value_threshold = 1e-130
		self.length_threshold = .85
		self.processors = None
		self.multi_hits = True

		self.ba = GOBlastFilter(self.my_blast_exe, self.my_blast_db, self.my_blast_file, self.e_value_threshold, self.length_threshold, blast_processors=self.processors, multi_hits = self.multi_hits)

		self.outfile_handle = open(self.my_blast_outfile)
		
		self.records = self.ba.runBlast(self.outfile_handle)

		self.filtered = self.ba.filterBlast(self.records)
		
	def testHitID(self):
		assert 124598 == self.filtered[str(3)][0].hit_id, 'wrong hit_id'
	def testQueryID(self):
		assert 3 == self.filtered[str(3)][0].query_id, 'wrong query_id'



if __name__ == "__main__":
            unittest.main()
    


