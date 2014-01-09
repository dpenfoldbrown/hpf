#!/usr/bin/python

import unittest
from hpf.hddb.db import Session, Family, Protein, Domain, PDBSeqRes
from hpf.pdb.mapper import ProteinToPDBAtom
from hpf.pdb.mapper import AlignmentToPDBMapper

class CompareToPDBMap(unittest.TestCase):

# Compare ProteinToPDBAtom mapper to PDBMap (make sure ProteinToPDBAtom is right)
# Compare on 2e32.

    def setUp(self, ):
        session = Session()
        self.fam  = session.query(Family).get(19187)
        self.prot = session.query(Protein).get(1171456)
        self.dom  = session.query(Domain).get(1307995)
        self.alignment = self.fam.alignment.alignment
        self.psr  = session.query(PDBSeqRes).filter_by(sequence_key=self.dom.parent_id[3:]).first()

        self.ptpdba = ProteinToPDBAtom(self.prot,self.dom,self.alignment, pdbseqres=self.psr, inverse=False)
        self.pdbmap = AlignmentToPDBMapper(self.fam.id, self.prot.id, self.dom.id, debug=True)

    def testMapIterate(self, ):
        print "Print ProteinToPDBAtom map"
        for col, atom in self.ptpdba:
            print col, atom

        print "Print AlignmentToPDBMapper map"
        for col in self.pdbmap.alignment_pdbatom_map.keys():
            print col, self.pdbmap.alignment_pdbatom_map[col]
        print "...Done"


if __name__ == "__main__":
    unittest.main()
