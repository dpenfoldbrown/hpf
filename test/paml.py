'''
Created on Sep 27, 2009

@author: patrick
'''
import os
import unittest
from hpf.paml import *
import test.hpf
import numpy

class PAMLTest(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def testCodeMLParser(self):
        """Test the parsing of ancestral sequences from an RST file"""
        input = os.path.join(os.path.dirname(test.hpf.__file__),"dat","codeml.rst")
        alignment = codeml_label(CodeMLParser().parse(input))
        assert len(list(alignment)) == 5
        assert [(record.id,str(record.seq)) for record in alignment] == [('seq1', 'ACC'), ('seq2', 'ACG'), ('seq3', 'AGC'), ('node4', 'ACC'), ('node5', 'ACC')]
        
        site_prob = {"seq1":numpy.array([1.0,1.0,1.0]),
                     "seq2":numpy.array([1.0,1.0,1.0]),
                     "seq3":numpy.array([1.0,1.0,1.0]),
                     "node4":numpy.array([0.978,0.478,0.999]),
                     "node5":numpy.array([1.0,1.0,1.0])
                     }
        seq_prob = {"seq1":1.0,
                    "seq2":1.0,
                    "seq3":1.0,
                    "node4":0.46740,
                    "node5":1.0
                    }

        #((seq2: 0.537243, seq1: 0.000004): 0.255741, seq3: 0.281503);
        branch_length = {"seq1":0.000004,
                         "seq2":0.537243,
                         "seq3":0.281503,
                         "node5":0.255741
                         }

        for branch in alignment.branches():
            desc = branch[1]
            assert branch.distance()==branch_length[desc.id]
        
        for record in alignment:
            assert all(record.letter_annotations["probability"] == site_prob[record.id])
            assert record.annotations["probability"] == seq_prob[record.id]

            

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
