#!/usr/bin/env python
import unittest
from hpf.scripts import domain_sccs
from hpf.pdb.scop import SCOPDomain

class test_domain_sccs(unittest.TestCase):

    def setUp(self):
        # Make a bunch of fake scop pdb entries
        self.domains = [SCOPDomain("1aa1","a.1.2",part_text=part) for part in ["A:,B:","D:118-211,D:372-450","E:1-384","E:385-422,E:431-487"]]
        for domain in self.domains:
            for part in domain.parts:
                part.chain_len = 30
        domain_sccs.pdb = {"1aa1":self.domains}

    def test_scop(self):
        """Test the scop parts and chain length stuffs"""
        self.failUnlessEqual(60,self.domains[0].length(), "Should sum the chains")
        self.failUnlessEqual((211-118+1)+(450-372+1),self.domains[1].length(), "Should sum the chains %i"%self.domains[1].length())
        
        

    def test_pdb_domain(self):
        domain_key = 0
        parent_key = 0
        domain_type = "fold_recognition"
        domain_len = 100
        
        task = domain_key, parent_key, domain_type, "1aa1", "A", 0, 30, domain_len
        result = domain_sccs.pdb_domain(task)
        self.failUnlessEqual(result[3],None, "Alignment overlap is NOT 50% and shouldn't count")

        task = domain_key, parent_key, domain_type, "1aa1", "A", 0, 31, domain_len
        result = domain_sccs.pdb_domain(task)
        self.failUnlessEqual(result[3],"a.1.2", "Alignment overlap IS 50%")

        task = domain_key, parent_key, domain_type, "1aa1", "E", 1, 487, domain_len
        result = domain_sccs.pdb_domain(task)
        self.failUnlessEqual(result[7],2, "Alignment overlap IS 30 AA, spans two domains!")

    def test_sort(self):
        import random
        psi_domain,mcm_domain = (0,1)
        parent_key = 2
        domain_types = ["psiblast","fold_recognition","msa","pfam","unassigned","user_defined"]

        psi_test = [(2,1),(1,0.8)]
        mcm_test = [(1,0.75),(1,0.8),(1,0.9)]
        psi_results = [(psi_domain, parent_key, random.sample(domain_types[1:2],1), "a.1.2", confidence, "1aa1","A", size) for size,confidence in psi_test]
        mcm_results = [(mcm_domain, parent_key, random.sample(domain_types[2:],1), "a.1.2", confidence, "1aa1","A", size) for size,confidence in mcm_test]
        
        # returns a dict of each domain_key to best result
        s = domain_sccs.sort(psi_results+mcm_results)
        self.failUnlessEqual(s[psi_domain][7], 2, "Two domain protein should be chosen")
        self.failUnlessEqual(s[mcm_domain][4], 0.9, "MCM top scrore should be chosen "+str(s[psi_domain][4]))
        


if __name__ == "__main__":
    unittest.main()
