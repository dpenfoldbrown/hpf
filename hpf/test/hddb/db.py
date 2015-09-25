'''
Created on Jan 6, 2010

@author: Patrick
'''
import unittest
from hpf.hddb.db import *

class Test(unittest.TestCase):


    def setUp(self):
        self.session = Session()
    
    def tearDown(self):
        self.session.close()

    def testName(self):
        print ""
        p = self.session.query(Protein).filter(Protein.experiment_key==880).limit(1).all()[0]
        print p
        print p.sequence
        print p.experiment
        print p.experiment.proteins


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()