#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
from hpf.pdb.cluster import BioPDBStruct
from hpf.hddb.db import Session, Family, Protein, Domain
from hpf.pdb.possel_cluster_driver import get_possel_sites, get_domain_struct, get_domain_sites


class PosSelClusterTest(unittest.TestCase):
    
    def setUp(self,):
        session = Session()
        #family_id  = 19187
        #protein_id = 1151960
        #domain_id  = 1302435
        family_id = 18883
        protein_id = 1063014
        domain_id = 1212855
        
        self.family = session.query(Family).get(family_id)
        self.protein = session.query(Protein).get(protein_id)
        self.domain  = session.query(Domain).get(domain_id)

        self.ps_sites = get_possel_sites(self.family)
        (self.pdb_id, self.pdb_chain, self.struct) = get_domain_struct(self.domain)
        
        print "Testing on Family {0}, Protein {1}, Domain {2}, Structure {3}, PDB {4}{5}".format(family_id, protein_id, domain_id, self.struct.id, self.pdb_id, self.pdb_chain)
        print "Sites of +Sel for family: ", self.ps_sites

        self.pdb_struct = BioPDBStruct(self.pdb_id, self.pdb_chain, debug=True)
        self.cluster_id_str = "Family {0}, Protein {1}, Domain {2}, Structure {3}".format(family_id, protein_id, domain_id, self.struct.id)


    def testRunClusterAnalysis(self, ):
        domain_sites = get_domain_sites(self.family, self.protein, self.domain, self.ps_sites)
        print "Local sites for domain {0}: ".format(self.domain.id), domain_sites

        self.pdb_struct.cluster_analysis(domain_sites, sample_size=10, store_file='ps-cluster-test.pkl', tag=self.cluster_id_str, report=True)


if __name__ == "__main__":
    unittest.main()
