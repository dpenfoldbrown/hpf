'''
Created on Nov 2, 2010

@author: Patrick
'''
import unittest
from hpf.mapper import MapperInterface

class TestMapper(MapperInterface):
    def __init__(self, *tuples, **kwargs):
        super(TestMapper,self).__init__(**kwargs)
        self.tuples = tuples
    def mapping(self):
        return iter(self.tuples)

class Test(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass

    def testInterface(self):
        m1 = TestMapper(('A','a'),('B','b'),('C','c'),('E','e'))
        m3 = TestMapper(('A','a'),('B','b'),('C','c'),('E','e'), inverse=True)
        m2 = TestMapper(('d','D'),('c','Z'),('a','A'))
        chain = MapperChain(m1,m2)
        chain._dict_()
        inv = chain.inverse()
        print "Forward",chain._forward
        print "Inverse",inv._forward
        assert chain['A']=='A'
        assert chain['C']=='Z'
        assert inv['Z']=='C'
        assert list(m1)!=list(m2)
        assert list(m1.inverse())==list(m3)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()