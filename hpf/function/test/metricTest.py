import unittest
from hpf.function.metric import *
from hpf.function.term import *

class MetricTest(unittest.TestCase):

	def setUp(self):
		mi = MutualInformation()
		mi.set_metric("F","A",0)
		mi.set_metric("F","B",1)
		mi.set_metric("F","C",2)
		mi.set_metric("F",'D',3)
		mi.set_metric("F","E",3)
		mi.set_metric("X","E",3)
		self.mi = mi

	def testGetAllMetric(self):
		met = self.mi.get_all_metric("X")
		assert 3 == met.get_metric("E",None)

	def testMI(self):
		terms = Terms()
		terms.append(Term("A"))
		terms.append(Term("B"))
		terms.append(Term("C"))
		max = terms.getFeaturedTerm(FuncTerm("F"), self.mi).get_id()
		assert max=="C", max

		terms = Terms()
		terms.append(Term("D"))
		terms.append(Term("E"))
		max = terms.getFeaturedTerm(FuncTerm("F"), self.mi).get_id()
		assert max=="D" or max=="E", max

		terms = Terms()
		max = terms.getFeaturedTerm(FuncTerm("F"), self.mi)
		assert max==None
		
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
