#!/usr/bin/python

import unittest
from hpf.pdb.mapper import PDBAtomSeqResMapper
from hpf.hddb.db import *

class Test2e32(unittest.TestCase):
    
    def setUp(self, ):
        self.session=Session()
        fam_id = 19187;
        prot_id = 1171456;
        dom_id  = 1307995;

        self.family = self.session.query(Family).get(fam_id)
        self.protein = self.session.query(Protein).get(prot_id)
        self.domain = self.session.query(Domain).get(dom_id)
        self.pdbseqres = self.session.query(PDBSeqRes).filter_by(sequence_key=self.domain.parent_id[3:], chain=self.domain.sccs.chain).first()
        self.pdbid = self.pdbseqres.pdb.pdbId+self.pdbseqres.chain

        self.map = PDBAtomSeqResMapper(self.pdbseqres, inverse=False)

    def testChain(self, ):
        print "Print whole map"
        for seq, atom in self.map:
            print seq, atom
        print "...Done"

if __name__ == "__main__":
    unittest.main()
