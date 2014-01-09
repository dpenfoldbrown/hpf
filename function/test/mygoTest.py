
import unittest
import os
import MySQLdb
from hpf.function.createdbs.mygo import Mygo
import hpf.function.term 

class SimpleTermTestCase(unittest.TestCase):
	def setUp(self):
		sql_user="kdrew"
		sql_password="kdrew_nyu"
		sql_host="localhost"
		sql_db = "mygo"

		self.mygo_conn = MySQLdb.connect(host=sql_host, user=sql_user, passwd=sql_password, db=sql_db)

	def testType(self):
		assert hpf.function.term.is_type("GO:0006139", "biological_process", self.mygo_conn)

class SimpleMygoRunTestCase(unittest.TestCase):
	def setUp(self):
		sql_user="kdrew"
		sql_password="kdrew_nyu"
		sql_host="mcpeepants.bio.nyu.edu"

		self.mg = Mygo(sql_user,sql_password, sql_host, e_codes =('TAS','IDA','IMP'), db_name="mygo")
		self.mg.connect()


	def testDBConnect(self):
		assert self.mg.isConnected(), 'did not connect to db'

	def testNumTerms(self):
		self.terms = self.mg.getGoTerms([1605014])
		self.terms.print_terms()
		print "number of terms: ", len(self.terms)
		assert 21 == len(self.terms), 'wrong number of terms'

	def testAllNumTerms(self):
		self.all_terms = self.mg.get_all_terms()
		assert 20608 == len(self.all_terms), 'wrong number of terms'


if __name__ == "__main__":
            unittest.main()
    


