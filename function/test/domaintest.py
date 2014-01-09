"""Unit test for domain.py
"""

import MySQLdb

import unittest
import domain.mcm_domain

class DomainTestCase(unittest.TestCase):
	def setUp(self):
		conn = MySQLdb.connect(host="mcpeepants.bio.nyu.edu", user="kdrew", passwd="kdrew_nyu", db="ecoli_benchmark")
		self.d3 = domain.mcm_domain.MCMDomain(37769,37769)
		self.d3.load_process(conn)
		self.d3.load_localization(conn)
		self.d3.load_function(conn)
		self.d3.load_mcm_structures(conn, scale=0.9)

	def testLengthsOfProcess(self):
		""" tests the lengths of process"""
		proc_len = len(self.d3.proc_terms)
		self.assertEqual(15, proc_len)
	def testLengthsOfStructure(self):
		""" tests the lengths of structure"""
		struct_len = len(self.d3.structures)
		self.assertEqual(6, struct_len)
	def testLengthsOfLocalization(self):
		""" tests the lengths of localization"""
		loc_len = len(self.d3.loc_terms)
		self.assertEqual(14, loc_len)
	def testLengthsOfFunction(self):
		""" tests the lengths of function"""
		func_len = len(self.d3.func_terms)
		self.assertEqual(3, func_len)



if __name__ == "__main__":
	unittest.main()
