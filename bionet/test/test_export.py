'''
Created on Sep 1, 2009

@author: patrick
'''
import unittest
import cStringIO
import random
import string
import networkx as nx
from hpf.bionet.export import XGMMLExporter
from hpf.bionet import Node,Edge

class XGMMLExporterTest(unittest.TestCase):
    """Doesn't really assert anything yet, just a sample use."""

    def setUp(self):
        self.output = cStringIO.StringIO()
        
    def tearDown(self):

        self.output.close()

    def testExport(self):
        graph = nx.MultiDiGraph(name="TestGraph")
        graph.graph["random"] = random.random()
        nodes = [Node("Node%i"%i) for i in range(10)]
        edges = [Edge((nodes[i],nodes[i+1]), {"random":random.random()}) for i in range(len(nodes)-1)]
        edges += [Edge((nodes[i+1],nodes[i]), {"random":random.random()}) for i in range(len(nodes)-1)]
        graph.add_nodes_from(nodes)
        for edge in edges:
            graph.add_edge(edge.get_source(),edge.get_target(),attr_dict=edge.attributes())
            #s,t = edge
            #print graph.get_edge_data(s,t)[0]
        print ""
        print "Graph",graph
        exporter = XGMMLExporter(self.output)
        exporter.write(graph)
        print self.output.getvalue()
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()