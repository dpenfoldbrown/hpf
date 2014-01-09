'''
Created on Mar 1, 2010

@author: patrick
'''
import unittest
from hpf.amnh.align import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Generic import Alignment
from Bio.Alphabet import IUPAC, Gapped

class Test(unittest.TestCase):


    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testCulledColumnMapper(self):
        align = Alignment(Gapped(IUPAC.protein, "-"))
        original = "ABCDEFGHI"
        align.add_sequence("test",original)
        culled = [0,1,4,8]
        # should yield
        result = "CDFGH"
        mapper = CulledColumnMapper(align,culled)
        for i,aa in enumerate(result):
            assert original[mapper[i]]==aa

    def testMap(self):
        alignment = SeedAlignment(Gapped(IUPAC.protein, "-"))
        seeds = [SeqRecord(Seq("A-BCD"), "seed1"),
                 SeqRecord(Seq("A-BCD"), "seed1")]
        target = SeqRecord(Seq("-ZBCD"), "target")
        alignment._records = seeds+[target]
        alignment._targets = [target]
        map = alignment.map(0)
        assert map==[None,1,2,3], map
                
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()