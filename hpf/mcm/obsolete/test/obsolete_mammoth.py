'''
Created on May 11, 2010

@author: patrick
'''
import unittest
import os
from hpf.hddb.mcm.mammoth import MammothMultiParser

DIR = os.path.dirname(__file__)

class TestParser(unittest.TestCase):


    def setUp(self):
        self.parser = MammothMultiParser()

    def tearDown(self):
        pass

    def testParser(self):
        file = os.path.join(DIR,"data.mammoth")
        with open(file) as handle:
            #result = MammothMultiParser.expression.parseFile(handle)
            result = self.parser.parse(handle)
        print result
        #d = result.asDict()
        #print d.values()
        
        scores = list(result)
        assert len(scores)==11
        print scores[0].__dict__
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()