"""Unit test for domain.py
"""

import MySQLdb

import unittest
import domain.mcm_domains

class DomainTestCase(unittest.TestCase):
	def setUp(self):
		conn = MySQLdb.connect(host="mcpeepants.bio.nyu.edu", user="kdrew", passwd="kdrew_nyu", db="ecoli_benchmark")
		self.ecoli_domains = domain.mcm_domains.MCMDomains("ecoli_benchmark", conn)

	def testLengthsOfDomains(self):
		""" tests the lengths of domain list"""
		domain_len = len(self.ecoli_domains)
		self.assertEqual(1352, domain_len)


if __name__ == "__main__":
	unittest.main()
