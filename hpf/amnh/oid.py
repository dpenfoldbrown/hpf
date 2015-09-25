'''
Created on Feb 24, 2010

@author: patrick
'''

from pyparsing import *
from hpf.parsing import *
import string
from string import punctuation
from hpf.runtime import runtime

from hpf.mapper import MapperInterface
class BioTreeNodeMapper(MapperInterface):
    """
    Map a biopython tree's indices to database keys.
    """
    def __init__(self,tree):
        """
        @param tree: Bio.Nexus.Trees.Tree generated from the database's
        TreeFactory.  Expects node's data object to contain an 'id' 
        (tree_node.id) attribute.
        """
        super(BioTreeNodeMapper,self).__init__()
        self.tree = tree

    def mapping(self):
        """
        @return generator yielding (biopython_node_id, tree_node_key)
        """
        # Index nodes by their database id's
        self.order_to_dbkey = {}
        for node_id in self.tree.all_ids():
            node = self.tree.node(node_id)
            tree_node_key = node.get_data().id
            yield (node_id,tree_node_key)

class DiagCharsParser(object):
    """
    Parse the diagnostic characters file.
    """
    def __init__(self, tree):
        """
        @param tree: The tree whose nodes will be used to store character info.
        @type tree: Bio.Nexus.Trees.Tree
        """
        self.tree = tree
        self.mapping = BioTreeNodeMapper(self.tree)
    
    def parse(self,handle):
        from hpf.hddb.db import TreeNodeDiagChars
        #tree = (poplar#1104386,poplar#1107391)
        #0 0:M,1:E,2:P,3:AG,4:IK
        #1 3:G,4:I,8:D,9:S,15:I,18:D,19:T,23:M
        #2 3:A,4:K,8:K,9:R,15:V
        while(True):
            line = handle.readline()
            if line.startswith("tree = "):
                break
        for line in handle:
            if line.strip() == "":
                continue
            try:
                node,chars = line.split()
            except:
                runtime().debug("Can't split line",line)
                raise
            tree_node_key = self.mapping[int(node)]
            for parts in chars.split(","):
                if parts.strip()=="":
                    continue
                try:
                    column,aa = parts.split(":")
                except:
                    runtime().debug("Can't split line",line)
                    runtime().debug("Can't split part",parts)
                    raise
                assert column.isdigit()
                column = int(column)
                # Indices here start at 1, map back to 0
                #print column, parts, line
                assert column >= 0
                assert len(aa) > 0
                yield TreeNodeDiagChars(
                    tree_node_key=tree_node_key,
                    column = column,
                    aa = aa)

class OIDTreeParser(object):
    """
    Parse the oid.tre PAUP run file for the given Nexus Tree.
    Uses the oid translation table to rename node taxon correctly.
    @deprecated
    """
    
    def __init__(self):
        assert False, "Dummy, use the NEXUS parser" 
    
    TRANSLATION = "Translate"
    NUMBER = "number"
    TAXON = "taxon"
    TREE = "tree"

    name_line = Group(Word(nums).setResultsName(NUMBER)+
                 Word(alphanums+punctuation).setParseAction(lambda tokens:[tokens[0].replace(",","")]).setResultsName(TAXON)
                 )
    
    tree = (Consume(Literal(";")).setResultsName(TREE)
            )
        
    translation = (OneOrMore(name_line).setResultsName(TRANSLATION)+
                   Literal(";")
                   )
    
    expression = (Consume(Literal(TRANSLATION)).suppress()+
                  translation+
                  Consume(Literal("tree PAUP_1 = [&R]")).suppress()+
                  tree+
                  empty
                  )
    
    
    def parse(self,handle):
        ps = OIDTreeParser
        result = OIDTreeParser.expression.parseFile(handle)
        from Bio.Nexus.Trees import Tree
        self.taxon = {}
        runtime().debug("Found",len(result[ps.TRANSLATION]),"species")
        print result
        for translation in result[ps.TRANSLATION]:
            print translation
            num = translation[ps.NUMBER]
            taxon = translation[ps.TAXON]
            runtime().debug(num,taxon)
            self.taxon[num] = taxon
        
        self.tree = Tree(result[ps.TREE])
        assert self.tree.count_terminals() == len(self.taxon.keys(), 
               "Translation count and # of tree terminals don't match")
        for node in self.tree.get_terminals():
            runtime().debug("Rename",node.get_data().taxon,"to",self.taxon[node.get_data().taxon])
            node.get_data().taxon = self.taxon[node.get_data().taxon]
        self.tree.no
        return self

def phylip(handle):
    seqs,columns = handle.readline().split()
    from Bio.Align.Generic import Alignment
    from Bio.Alphabet import IUPAC, Gapped
    alignment = Alignment(Gapped(IUPAC.protein, "-"))
    for line in handle:
        name,seq = line.split()
        alignment.add_sequence(name, seq)
    return alignment
        
def oid_amnh():
    """@deprecated"""
    from hpf.amnh import DATA_FOLDER
    import os
    from collections import defaultdict
    file = os.path.join(DATA_FOLDER,"oid.amnh")
    map = defaultdict(lambda: None)
    with open(file) as handle:
        for line in handle:
            oid,amnh = line.split()
            map[oid] = amnh
    return map

def oid_hpf():
#    from hpf.amnh import DATA_FOLDER
#    import os
#    from collections import defaultdict
#    file = os.path.join(DATA_FOLDER,"oid.hpf")
#    map = defaultdict(lambda: None)
#    with open(file) as handle:
#        for line in handle:
#            oid,hpf = line.split()
#            map[oid] = hpf
#    return map
    class lookup(dict):
        def __init__(self):
            self._oid_amnh = oid_amnh()
            
        def __getitem__(self, key):
            from hpf.hddb.db import Session, SequenceAc
            amnh = self._oid_amnh[key]
            session = Session()
            try:
                ac = session.query(SequenceAc).filter(SequenceAc.ac==amnh).first()
                return ac.protein_key if ac else None
            finally:
                session.close()
            
    return lookup()


def index(records):
    """
    In-place rename's record id's from OID name to HPF sequence_key
    """
    from hpf.hddb.db import Session, SequenceAc
    mapping = oid_amnh()
    session = Session()
    for record in records:
        amnh = mapping[record.id]
        id = str(session.query(SequenceAc).filter(SequenceAc.ac==amnh).one().sequence_key)
        record.id = id
     
