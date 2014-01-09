import os
DATA_FOLDER = os.path.join(os.path.dirname(__file__),"data")
JNLP = os.path.join(DATA_FOLDER,"default.jnlp")

class Attributes():
    """
    Abstract Static attributes
    """
    VALUE = "Value"
    NAME = "Name"
    ACCESSION = "Accession"
    PLS_LLR = "pls_llr"
    BASE_LLR = "base_llr"
    DOMAINS = "Domains"
    COVERAGE = "Coverage"
    PDB = "PDB"
    ROSETTA = "Rosetta"
    QUALITY = "Quality"
    SUPERFAMILIES = "Superfamilies"
    GINZU = "Ginzu"
    HDDB = "hddb"
    URL = "Human Proteome Folding Project URL"
    HDDB_URL = "HPF Report URL"
    
    IEA_ACC = "IEA acc"
    IEA_NAME = "IEA name"
    IEA_LLR = "IEA LLR(acc)"
    PMID = "Pubmed Citations"

    TYPES = {int:"integer",
             basestring:"string",
             str:"string",
             float:"real",
             long:"integer",
             bool:"boolean"
             }

    ROW_TO_COLUMN = {"size":DOMAINS,
                     "pdb":PDB,
                     "rosetta":ROSETTA,
                     "quality":QUALITY,
                     "sccs":SUPERFAMILIES,
                     "ginzu":GINZU,
                     "mf_value":VALUE,
                     "mf_name":NAME,
                     "mf_acc":ACCESSION,
                     "mf_pls_llr":PLS_LLR,
                     "mf_base_llr":BASE_LLR,
                     "coverage":COVERAGE,
                     "parent_sequence_key":HDDB,
                     "url":URL,
                     "hddb_url":HDDB_URL
                     #"iea_acc":IEA_ACC,
                     #"iea_name":IEA_NAME,
                     #"iea_llr":IEA_LLR,
                     #"pmid_count":PMID
                     }
    
    HPF_FUNCTION_INFO = "HPF Function";
    HPF_STRUCTURE_INFO = "HPF Structure";
    PREFIX = {
              VALUE: HPF_FUNCTION_INFO,
              NAME: HPF_FUNCTION_INFO,
              ACCESSION: HPF_FUNCTION_INFO,
              PLS_LLR: HPF_FUNCTION_INFO,
              BASE_LLR: HPF_FUNCTION_INFO,
              DOMAINS: HPF_STRUCTURE_INFO,
              COVERAGE: HPF_STRUCTURE_INFO,
              PDB: HPF_STRUCTURE_INFO,
              ROSETTA: HPF_STRUCTURE_INFO,
              QUALITY: HPF_STRUCTURE_INFO,
              SUPERFAMILIES: HPF_STRUCTURE_INFO,
              GINZU: HPF_STRUCTURE_INFO}

    #def __init__(self,attributes={}):
    #    self._dict = defaultdict(lambda: None)
    #    for col in attributes:
    #        self.set_attribute(col,attributes[col])
    
    def __init__(self,attr=None):
        self._dict = dict()
        for key in attr:
            self.set(key, attr[key])

    def keys(self):
        return self._dict.keys()
    
    def values(self):
        return self._dict.values()
    
    def has_key(self,key):
        return self._dict.has_key(key)
    
    def get(self, key):
        return self._dict[key]
    
    def set(self,key,value):
        """
        Set a node attribute, converting the column to HPF columns if necessary.
        """
        if Attributes.ROW_TO_COLUMN.has_key(key):
            key = Attributes.ROW_TO_COLUMN[key]
            if Attributes.PREFIX.has_key(key):
                key = Attributes.PREFIX[key]+": "+key
        self._dict[key] = value
    
    @staticmethod
    def get_type(value):
        return Attributes.TYPES[type(value)]
    
    def attributes(self):
        return self._dict
        
class Node(Attributes):
    def __init__(self,name,attr=None):
        Attributes.__init__(self,attr)
        self.set("name",str(name))
    
    def label(self):
        return self.get("name")
    
    #def __hash__(self):
    #    return hash(self["name"])

class Edge(Attributes,tuple):

    def __new__(cls, (source, target) ,attr=None):
        return tuple.__new__(cls, (source, target))

    def __init__(self,(source,target),attr=None):
        tuple.__init__((source, target))
        Attributes.__init__(self,attr)
        self.source = source
        self.target = target
        
    def label(self):
        inter_type = " ("+self.get("interactionType")+") " if self.has_key("interactionType") else " "
        return self.get_source().label()+inter_type+self.get_target().label()
    #def __hash__(self):
    #    hash((self.source,self.target))
    
    def get_source(self):
        return self.source

    def get_target(self):
        return self.target
