'''
Created on Nov 3, 2010

@author: Patrick
'''
import unittest
from pyparsing import *
from hpf.parsing import Consume
from hpf.amnh.codeml import PositiveSelectionParser
from hpf.test.amnh import DIR
import os

LRT_FILE = os.path.join(DIR,"codeml.lrt") 

class Test(unittest.TestCase):


    def setUp(self):
        self.file = LRT_FILE 

    def tearDown(self):
        pass


    def testParserSmallFile(self):
        parser = PositiveSelectionParser()
        with open(self.file+"1") as handle:
            models = list(parser.parse(handle))
        assert len(models)==1
        assert len(models[0].positive_selection)==4

    def testParserBigFile(self):
#        ps = PositiveSelectionParser
#        exp = Consume(ps.selection_header).suppress()+OneOrMore(ps.selection_line)+empty
#        exp.setDebug(True)
#        with open(self.file+"1") as handle:
#            r = exp.parseFile(handle)
#        print "result",r
        parser = PositiveSelectionParser()
        with open(self.file) as handle:
            models = list(parser.parse(handle))
        assert len(models)==3
        
        index = {}
        for model in models:
            index[model.model] = sorted([(ps.column, ps.probability) for ps in model.positive_selection])
            
        assert len(index[2]) == 2
        assert index[2] == [(188,float("0.573")),
                            (193,float("0.560"))] 

        assert len(index[8]) == 3 
        assert index[8][-1] == (582,float("0.932"))

        print index[3]
        assert len(index[3]) == 4
        assert index[3] == [(2,float("0.925")),
                            (3,float("0.998")),
                            (4,float("0.917")),
                            (582,float("1.0"))] 


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()