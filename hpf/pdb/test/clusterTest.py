
import os
import MySQLdb
import unittest

import hpf.pdb.cluster as cluster

class SimpleClusterTestCase(unittest.TestCase):

	def setUp(self):
		#kdrew: these cases are made up
		#kdrew: two tight clusters
		self.distance_list1 = [1,1,1,1,1,1,9,10,11,9,10,11,9,10,11]
		#kdrew: one tight cluster
		self.distance_list2 = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
		#kdrew: no cluster
		self.distance_list3 = [20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]

	def testInverseDistanceMean(self):
		print "testing inverse_distance_mean"
		Ms1 = cluster.inverse_distance_mean(self.distance_list1)
		Ms2 = cluster.inverse_distance_mean(self.distance_list2)
		Ms3 = cluster.inverse_distance_mean(self.distance_list3)
		assert Ms1 > Ms3, "no-cluster sample has larger Ms value than two-clustered sample"

		assert Ms2 == 0.5, "Ms of clustered sample is not correct"


class StructureClusterTestCase(unittest.TestCase):

	def setUp(self):
		#kdrew: load pdb file
		pdb_filename = "data/1UBQ.pdb"
		pdb_handle = open(pdb_filename)
		pdb_string = pdb_handle.read()

		#kdrew: using list of residue ids that are near each other
		clustered_ids =  [2,3,15,16,64]
		self.dist_list = cluster.structure_distance_calc(pdb_string, clustered_ids)
		self.Ms = cluster.inverse_distance_mean(self.dist_list)
		self.random_Ms_list = cluster.random_structure_distance_calc(pdb_string, len(clustered_ids))
		print self.random_Ms_list

	def testInverseDistanceMean(self):
 		distance_list = map(round, [5.8570064879595272, 7.0676930465322281, 5.8547376542420757, 6.1065838240377914, 5.2907485292725829, 8.1445831078085273, 6.9240371893859729, 5.6990176346454655, 11.13853720198483, 11.6613908690173],[4]*10)

		print distance_list
		print self.dist_list
		assert map(round,self.dist_list,[4]*len(self.dist_list)) == distance_list, "distance list does not match"
		assert round(self.Ms,6) == round(0.145399896701,6) , "Ms does not match"

if __name__ == "__main__":
            unittest.main()
    


