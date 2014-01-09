
import unittest
import MySQLdb

from hpf.function.createdbs.astral import Astral, SuperfamilyEntry


class SimpleAstralRunTestCase(unittest.TestCase):
	def setUp(self):
		sql_user="kdrew"
		sql_password="kdrew_nyu"
		#sql_host="mcpeepants.bio.nyu.edu"
		sql_host="127.0.0.1"
		db_name = "pdb"
		astral_table = "astral95_1_75"

		self.astl = Astral(sql_user,sql_password, sql_host, db_name=db_name, table=astral_table)
		self.astl.connect()

		self.sfam = SuperfamilyEntry(sf_id = 'a.1.1')



	def testDBConnect(self):
		assert self.astl.isConnected(), 'did not connect to db'

	def testAllNumSFs(self):
		self.all_superfamilies= self.astl.get_all_superfamilies()
		assert 1961 == len(self.all_superfamilies), 'wrong number of terms'
	def testSF(self):
		assert "a.1.1" == self.sfam.get_id(), 'wrong superfamily id'



if __name__ == "__main__":
            unittest.main()
    


