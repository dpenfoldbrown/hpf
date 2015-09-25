'''
Created on Nov 3, 2010

@author: Patrick
'''
import unittest
from StringIO import StringIO
from Bio.Nexus.Trees import Tree
from hpf.amnh.oid import DiagCharsParser
from collections import defaultdict
diag_str = """
tree = (poplar#1104386,poplar#1107391)
0 0:M,1:E,2:P,3:AG,4:IK
1 3:G,4:I,8:D,9:S,15:I,18:D,19:T,23:M
2 3:A,4:K,8:K,9:R,15:V
"""
tree_str = "(poplar#1104386,poplar#1107391)"
class Test(unittest.TestCase):


    def setUp(self):
        self.io = StringIO(diag_str)
        self.io.seek(0)
        self.tree = Tree(tree_str)
        self.tree.node(0).get_data().id = 10
        self.tree.node(1).get_data().id = 20
        self.tree.node(2).get_data().id = 30
        
    def tearDown(self):
        pass

    def testParser(self):
        parser = DiagCharsParser(self.tree)
        diagchars = list(parser.parse(self.io))
        index = defaultdict(lambda: list())
        for d in diagchars:
            index[d.tree_node_key].append(d)
            index[d.tree_node_key].sort(cmp=lambda x,y:cmp(x.column,y.column))

        assert len(index[10]) == 5
        assert len(index[20]) == 8
        assert len(index[30]) == 5
        assert index[10][0].column == 0
        assert index[10][4].aa == "IK"
        assert index[20][2].column == 8 
        assert index[20][2].aa == "D"

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()