'''
Created on Aug 31, 2009

Export graphs for BioNetBuilder.

@author: patrick
'''
from xml import sax
from xml.sax import saxutils
from hpf.bionet import Attributes, Edge, Node
from collections import defaultdict
import networkx as nx

class NetworkXReader(sax.xmlreader.XMLReader):
    """
    Reads a NetworkX graph to produce SAX events.
    """
    
    def __init__(self):
        sax.xmlreader.XMLReader.__init__(self)
        self.depth = 0
    
    def _startElement(self,name,attrs=None):
        handler = self.getContentHandler()
        # For pretty printing
        handler.ignorableWhitespace(" "*self.depth)
        handler.startElement(name,attrs)
        self.depth += 1
        handler.ignorableWhitespace('\n')
        
        
    def _endElement(self,name):
        handler = self.getContentHandler()
        # For pretty printing
        self.depth -= 1
        handler.ignorableWhitespace(" "*self.depth)
        handler.endElement(name)
        handler.ignorableWhitespace('\n')

    
    def parse(self, source):
        """Parse a graph into XML events."""
        if not isinstance(source, NetworkXInputSource):
            raise NetworkXInputSourceException("Not a NetworkXInputSource")
        
        self._source = source
        handler = self.getContentHandler()
        handler.startDocument()
        self._graph(self._source.graph())        
        handler.endDocument()

    def _graph(self, graph):
        attrs = {"label":graph.name,
                 "xmlns:dc":"http://purl.org/dc/elements/1.1/",
                 "xmlns:xlink":"http://www.w3.org/1999/xlink",
                 "xmlns:rdf":"http://www.w3.org/1999/02/22-rdf-syntax-ns#",
                 "xmlns:cy":"http://www.cytoscape.org",
                 "xmlns":"http://www.cs.rpi.edu/XGMML",
                 "directed":"1"}
        self._startElement("graph",attrs)
        # Graph data... weird it's named 'graph'
        graph_data = graph.graph
        for attr in graph_data:
            value = graph_data[attr]
            self._att(attr, value)
        # Graph nodes
        for node in graph.nodes_iter():
            self._node(node)
        # Graph edges. nx blows away class, these are tuples
        for edge in graph.edges_iter():
            self._edge(edge)
        self._endElement("graph")

    def _att(self, name, value):
        #<att type="string" name="vizmap:Human Mystery NODE_SHAPE" value="ELLIPSE"/>
        if value==None:
            return
        try:
            attrs = {"type":Attributes.get_type(value),
                     "name":str(name),
                     "value":str(value)
                     }
        except KeyError:
            print name,value
            raise
        self._startElement("att", attrs)
        self._endElement("att")

    def _node(self, node):
        """Process a single node"""
        #<node label="21327701" id="-88">
        node_id = str(self._source.get_id(node))
        node_attrs = {"label":node.label(),
                      "id":node_id}
        self._startElement("node", node_attrs)
        for attr in node.keys():
            value = node.get(attr)
            self._att(attr, value)
        self._endElement("node")
        
    def _edge(self, edge):
        """Process a single edge"""
        #<edge label="GI:4506787 (PP.0.1) 4757960" source="-13" target="-82">
        source, target = edge
        # Odd hack.  I don't know why the data is keyed on 0, maybe a 
        # directional thing?
        attrs = self._source.graph().get_edge_data(source,target)[0]
        edge = Edge(edge,attrs)
        source_id = str(self._source.get_id(source))
        target_id = str(self._source.get_id(target))
        edge_attrs = {"label":edge.label(),
                      "source":source_id,
                      "target":target_id}
        self._startElement("edge", edge_attrs)
        for attr in edge.keys():
            value = edge.get(attr)
            self._att(attr, value)
        self._endElement("edge")


class XGMMLExporter():
    
    def __init__(self, file):
        """
        @param file: The file to write the XML to.
        @type file: Either string or open file_handle.
        """
        self._close = isinstance(file, basestring)
        self._file = file
    
    def __enter__(self):
        if self._close:
            self._file = open(self._file,"w")
        return self
        
    def __exit__(self, type, value, traceback):
        if self._close:
            self._file.close()
        return False
    
    def write(self, graph, pretty_print=True):
        """
        Write a graph to XML.
        @type graph: NetworkX graph type. 
        """
        # Use a graph/xgmml parser
        reader = NetworkXReader()
        source = NetworkXInputSource(graph)
        handler = XMLGenerator(self._file,pretty_print)
        reader.setContentHandler(handler)
        reader.parse(source)
    
    def _node(self,node):
        for attribute in node:
            name = XGMMLExporter.__attribute_name(attribute)
            value = XGMMLExporter.__atrribute_value(node, attribute)
            
    def _edge(self):
        pass        
    
            
    def __attribute_name(attribute):
        return Attributes.PREFIX[attribute] + ": " + attribute
    
    def __atrribute_value(node, attribute):
        return str(node.get_attribute(attribute))

class XMLGenerator(saxutils.XMLGenerator):

    def __init__(self, file, pretty_print=True):
        """
        @param pretty_print: Whether to print white-space and line breaks.
        """
        saxutils.XMLGenerator.__init__(self,file)
        self._pretty = pretty_print

    def ignorableWhitespace(self,s):
        if self._pretty:
            self.ignorableWhitespace(s)

def create_parser():
    return NetworkXReader()

class NetworkXInputSource(sax.xmlreader.InputSource):
    """
    Implements a NetworkX graph as an InputSource for the SAX interfaces.
    """
    def __init__(self, graph):
        if not isinstance(graph, nx.Graph):
            raise NetworkXInputSourceException("Not a NetworkX Graph")
        self._graph = graph
        self._unique_id = dict()
        self._unique_class = defaultdict(int)
    
    def graph(self):
        return self._graph
    
    def get_id(self,element):
        """
        Get a unique id for each type.
        """
        if not self._unique_id.has_key(element):
            # Count the number per class, assign next id to element
            self._unique_class[element.__class__] += 1
            self._unique_id[element] = self._unique_class[element.__class__]
        return self._unique_id[element]
        
class NetworkXInputSourceException(sax.SAXException):
    pass

if __name__=="__main__":
    pass
            
