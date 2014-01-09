'''
Created on Nov 11, 2009

@author: patrick
'''
import unittest
import os
from test.hpf import DATA_FOLDER
from hpf.pdb.psipred import parse

class TestPsipred(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def testParse(self):
        with open(os.path.join(DATA_FOLDER,"polypeptide.horiz")) as handle:
            prediction = parse(handle)
        
        assert str(prediction)=="CCCCCCCCCCCHHHHHHHHHHCCCCHHHHHHHHHCCCCCCCHHHHHHHHHHHHHHHCHHHHHCCCEEEEEECC", str(prediction)
        assert prediction.weights == [int(c) for c in "9987677677899999999408984999999860256888889999999999997293672317856887339"], prediction.weights
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()