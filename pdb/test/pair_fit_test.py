#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import hpf.pdb.pairfit as pf

class PairFitTest(unittest.TestCase):
    
    def setUp(self,):
		apdb_filename = "./data/p53_helix.pdb"
		structA_vector_ids = [ 19, 23, 26 ]

		bpdb_filename = "./data/oop_dimer_llll.pdb"
		structB_vector_ids = [ 1, 2, 3, 4 ]
		selection_id_type = "RESIDUE"
	
		self.pfit = pf.PairFit( apdb_filename, bpdb_filename, structA_vector_ids, structB_vector_ids, selection_id_type )	

    def testPairFit(self, ):
		results = self.pfit.pair_fit()
		print results[0]
		assert( results[0][2], 1.9310051202774048 ), "best match rmsd does not match expected"

if __name__ == "__main__":
    unittest.main()
