'''
Created on May 13, 2010

@author: Patrick
'''
from hpf.utilities import Struct
from hpf.runtime import runtime

class StructParser(object):
    """
    Simple domlette parsing into Struct's for nodes/values to a single depth.
    """
    
    @staticmethod
    def parse(dom, xpath, type=Struct):
        """
        Simple parsing of xpath to a depth of one to produce one of type.
        """
        e = dom.xpath("/mcm/param")[0]
        return StructParser.parseElement(e, type)

    @staticmethod
    def parseElement(element, type=Struct):
        """
        With a domlette element, parse into type.
        """
        from Ft.Xml.cDomlette import Element
        parameters = {}
        for node in element:
            if not isinstance(node, Element) :
                # skipping text nodes
                continue
            param_name = str(node.tagName.strip())
            #runtime().debug(param_name)            
            param_value = node.firstChild.nodeValue if node.firstChild != None else ''
            parameters[param_name] = param_value
            
        return type(**parameters)

class ClusterCenter(object):
    
    @staticmethod
    def parse(dom, xpath="/mcm/extract"):
        """
        Parse a DOM for all Cluster Center objects.
        @param xpath: default '/mcm/extract'
        @return: dictionary of ClusterCenter objects keyed by index.
        """
        from Ft.Xml.cDomlette import Element
        e = dom.xpath(xpath)[0]
        centers = {}
        import re
        center_name = re.compile("center[0-9]{3}$")
        center_size = re.compile("center[0-9]{3}size$")
        
        rank=1
        for node in e:
            if not isinstance(node, Element) :
                # skipping text nodes
                continue
            if not node.tagName.startswith("center"):
                # just avoid problems
                continue

            name = node.tagName.strip()
            value = node.firstChild.nodeValue
            index = int(value)
            # Care should be taken with these "num=" lines.
            # This assumes that they all read center001... ie "center%3i"
            if center_name.match(name):
                num = int(name[6:])
                if not centers.has_key(num):
                    centers[num] = ClusterCenter(index=num,dom=dom)
                centers[num].index = int(value)
            elif center_size.match(name):
                num = int(name[6:-4])
                if not centers.has_key(num):
                    centers[num] = ClusterCenter(index=num,dom=dom)
                centers[num].size = int(value)
                centers[num].rank = rank
                rank += 1
        return centers
    
    def __init__(self, 
                 index=None, 
                 size=None,
                 rank=None,
                 dom=None):
        self.index = index
        self.size = size
        self.rank = rank
        self.dom = dom
        
    def __iter__(self):
        """
        @return: generator of mcm entries for this cluster.
        """
        entries = self.dom.xpath("/mcm/confidence/data/entry[cluster_center_index=%i]" % self.index)
        from itertools import imap
        from hpf.hddb.db import McmData
        return imap(lambda e: StructParser.parseElement(e,type=McmData),entries)
        
    def __str__(self):
        return "<ClusterCenter:%i size:%i>" % (self.index, self.size)
    
    def __cmp__(self, other):
        return cmp(self.index, other.index)
    
    def entries(self):
        """
        @return: a list of mcm entries for this cluster center.
        @see: ClusterCenter.__iter__
        """
        return list(self)
    
    def decoy(self):
        """
        Get this cluster center's decoy atom record from the domlette.
        @return: str
        """
        #<decoys><decoy><name>decoy_10043.pdb</name>
        decoy_name = "decoy_%i.pdb" % self.index
        node = self.dom.xpath("/mcm/decoys/decoy[name=\"%s\"]" % decoy_name)[0]
        record = node.xpath(".//atomrecord")[0].firstChild.nodeValue
        return record
            
class McmLogFile(object):
    """MCM result log file, parses cluster centers into McmData parent type."""
    
    def __init__(self, dom):
        """Parse a logfile into McmData format"""
        self.dom = dom
        
    def __iter__(self):
        """Return the logfiles entries by parsing the file."""
        return self.entries()
    
    @staticmethod
    def parseString(st, **kwargs):
        from Ft.Xml.Domlette import NonvalidatingReader
        st = st.replace(">&", "&gt;&amp;")
        st = st.replace("<MAMMOTH>", "&gt;MAMMOTH&lt;")
        dom = NonvalidatingReader.parseString(st, **kwargs)
        return McmLogFile(dom)
        
    @staticmethod
    def parseFile(file):
        """File can be an open handle or filesystem path"""
        from Ft.Xml.Domlette import NonvalidatingReader
        if isinstance(file, basestring):
            dom = NonvalidatingReader.parseUri("file:%s" % self._file)
        else:
            dom = NonvalidatingReader.parseStream(file, **kwargs)
        return McmLogFile(dom)

    def convergence(self):
        conv = self.dom.xpath("/mcm/convergence")[0]
        return StructParser.parseElement(conv)

    def entries(self):
        """Returns McmData objs for each entry in the xml"""
        entries = self.dom.xpath("/mcm/confidence/data/entry")
        from itertools import imap
        from hpf.hddb.db import McmData
        return imap(lambda e: StructParser.parseElement(e,type=McmData),entries)

    def params(self):
        params = self.dom.xpath("/mcm/param")[0]
        return StructParser.parseElement(params)
    
    def cluster_centers(self):
        """
        @return: list of ClusterCenter objects extracted from the DOM.
        """
        return ClusterCenter.parse(self.dom).values()
        
