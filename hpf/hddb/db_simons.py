'''
Created on Jan 5, 2010

@author: patrick

Updated Apr 22, 2016 (jesus), to fit with HPF database imported
to Simons Foundation
'''
import sqlalchemy.ext.declarative
from sqlalchemy import create_engine, ForeignKeyConstraint
from sqlalchemy import *
from sqlalchemy.types import *
from sqlalchemy.orm import *
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.sql import and_, or_, select, func
from sqlalchemy.orm.interfaces import *
from hpf.runtime import runtime
from sqlalchemy.sql.expression import desc
from sqlalchemy.exc import IntegrityError
from datetime import datetime

# url = "mysql://dpb:dpb_nyu@ms2.bio.nyu.edu:3306"
# url = "mysql://dpb:dpb_nyu@handbanana.bio.nyu.edu:3306/"
# url = "mysql://kdrew:bonneau@127.0.0.1:13307/"
# url = "mysql://pfp:pfp_nyu@db/"
url = "mysql://nygo:189213955@db0.scdanet.org:3306/"

temp = None
engine = None
Base = None
Session = None
def clear():
    global engine, Session, Base
    runtime().debug("Clearing HDDB engine/session mapping")
    engine, Session, Base = (None,None,None)

def setup(db='hpf_original'):
    # create_engine params: pool_recycle=300, echo=True|"debug", pool_timeout=10)
    e = create_engine(url+db, echo=False, pool_recycle=1800)
    b = sqlalchemy.ext.declarative.declarative_base(bind=e)
    s = sessionmaker(bind=e)
    return e,b,s

def setup_scoped(db='hpf_original'):
    """
    Description of scoped session:
    http://docs.sqlalchemy.org/en/rel_0_7/orm/session.html#unitofwork-contextual
    """
    global engine, Base
    if engine == None:
        engine = create_engine(url+db, echo=False, pool_recycle=1800)
    if Base == None:
        Base = sqlalchemy.ext.declarative.declarative_base(bind=engine)
    s = scoped_session(sessionmaker(bind=engine))
    return s 

def rebind(**kwargs):
    global engine, Base, Session
    clear()
    engine, Base, Session = setup(**kwargs)

from hpf.hddb import tunnel
engine, Base, Session = setup()
ScopedSession = setup_scoped()


def push_to_db(session, object, exception_str=None, raise_on_duplicate=True):
    """
    A simple utility function to do the session add-flush-catch-refresh blocks.
    If raise_on_duplicate is true, will raise exception on IntegrityError. If false,
    will rollback and continue
        session - the Session to push ORM objects into
        object  - the ORM object to add to the session (and DB)
        exception_str   - a string to print when flushing the session fails (DB error)
    """
    if not exception_str:
        exception_str = "Failed to add object {0} to DB".format(object)

    session.add(object)
    try:
        session.flush()
    except IntegrityError:
        print "Given object already exists in DB"
        if raise_on_duplicate:
            session.rollback()
            print exception_str
            raise
        else:
            print "Rolling back, returning None"
            session.rollback()
            return None
    except Exception:
        session.rollback()
        print exception_str
        raise
    #finally:
    #    session.close()
    
    session.refresh(object)
    return object


from hpf.utilities import Struct
class BaseStruct(object):
    
#    _descriptor = {}
    
    def __init__(self, *args, **kwargs):
        #Base.__init__(self)
        self._update_vars(kwargs)
        
    def _update_vars(self, kwargs):
        for key in kwargs:
            self.__setattr__(key,kwargs[key])
            
    def __repr__(self):
        vars = [str(key)+":"+str(self.__dict__[key]) for key in self.__dict__ if not key.startswith("_")]
        return "<%s %s>" % (self.__class__.__name__," ".join(vars))

    def _set_text(self, text):
        self.compressed = func.compress(text)
    def _get_text(self):
        return column_property(select(["uncompress(compressed)"]),deferred=True)
    text = property(_get_text, _set_text)


#    def __setattr__(self,instance,value):
#        """
#        Executes a function for specific attributes, as defined by the Class'
#        _descriptor dictionary.  Useful for setting compressed fields from
#        text attributes.
#        """
#        Base.__setattr__(self,instance,value)
#        if hasattr(self.__class__, "_descriptor") and self.__class__._descriptor.has_key(instance):
#            self.__class__._descriptor[instance](self,instance,value)
            
class MammothFactory(object):
    
    def create(self, mammoth_score):
        ms = mammoth_score
        return Mammoth(zscore=ms.zscore,
                       ini_psi=ms.psi1,
                       end_psi=ms.psi2,
                       score=ms.score,
                       evalue=ms.evalue,
                       ln_e=ms.ln_e,
                       num_p=ms.num_p,
                       num_e=ms.num_e,
                       nss=ms.nss,
                       nsup=ms.nsup,
                       p_contact_order=ms.prediction_contact_order,
                       e_contact_order=ms.experiment_contact_order
                       )
    
class TreeNodeFactory(object):
    """
    Generate TreeNodes from a Bio.Nexus.Tree
    Ancestor is reference to ancestor TreeNode object.
    """
    def create(self, tree):
        from collections import defaultdict
        nodes = defaultdict(lambda: None)
        for id,node in tree.chain.items():
            tree_node = TreeNode()
            nodes[id] = tree_node
            
            data = node.get_data()
            tree_node.ancestor = nodes[node.get_prev()]
            tree_node.protein_key = data.taxon
            tree_node.branchlength = data.branchlength
            tree_node.branchlength_sum = tree.sum_branchlength(node=id)
            tree_node.tree_key = tree.name
            yield tree_node

class TreeFactory(object):
    """
    Create a Bio.Nexus.Tree object from a set of TreeNodes.
    Attaches 'diagchars' list and 'id' (tree_node.id) to node's data.
    """
    def __init__(self,
                 name_func=lambda node:str(node.id), ):
        self._name_func = name_func

    def create(self, nodes, tree_key=None):
        from Bio.Nexus.Trees import Tree,NodeData
        from Bio.Nexus.Nodes import Node
        map = {}
        root = None
        # Hopefully sorted ensures all of the nodes predecessors exist when building
        sorted_nodes = sorted(nodes, cmp=lambda x,y: cmp(x.id,y.id))
        tree = Tree(name = tree_key if tree_key else sorted_nodes[0].tree_key)
        
        # Create NodeData objects for each TreeNode object and index
        print len(sorted_nodes)
        for node in sorted_nodes:
            if node.ancestor_node==None:
                tree_node = tree.node(tree.root)
                data = tree_node.get_data()
                data.id = node.id
                data.diagchars = node.diagchars
            else:
                data = NodeData(taxon=self._name_func(node),
                     branchlength=node.branchlength,
                     support=None, 
                     comment=None)
                data.diagchars = node.diagchars
                data.id = node.id
                tree_node = Node(data)
                tree.add(tree_node,prev=map[node.ancestor_node].get_id())
            map[node.id] = tree_node
        return tree
        
class AlignmentFactory(object):
    """
    Create a Bio.Align.Generic object from a set of family_proteins.
    @param name_func: A function that will return a name for a record to be used in the alignment.
    """
    def __init__(self,
                 name_func=lambda p:str(p.protein_key)):
        self._name_func = name_func

    def create(self, *proteins, **kwargs):
        from Bio.Align.Generic import Alignment
        from Bio.Alphabet import IUPAC, Gapped
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        alignment = Alignment(Gapped(IUPAC.protein, "-"), **kwargs)
        for protein in proteins:
            record = SeqRecord(Seq(protein.sequence),id=self._name_func(protein),description=str(protein.protein_key))
            alignment._records.append(record)
        return alignment

class SeqRecordFactory(object):
    """Create a Bio.SeqRecord from an hddb Sequence object"""
    def create(self, *sequences, **kwargs):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        for sequence in sequences: 
            yield SeqRecord(Seq(sequence.sequence),str(sequence.id), **kwargs)
 
class SequenceFactory(object):
    """Create an hph.hddb.db.Sequence object from a SeqRecord, using record.id if is int()"""
    def create(self, record):
        from hashlib import sha1
        sequence = Sequence()
        sequence.sequence = str(record.seq)
        sequence.sha1 = sha1(str(record.seq)).hexdigest()
        try:
            sequence.id = int(record.id)
        except:
            pass
        return sequence
 
class Sequence(BaseStruct, Base):
    __tablename__ = 'sequence'
    __table_args__ = (
            ForeignKeyConstraint(['id'], ['sequenceAc.sequence_key']),
            ForeignKeyConstraint(['id'], ['hddb_iea_golite_062009.sequence_key']),
            {'autoload':True}
            )
    
    common_name = {'sequence':'Protein Sequence'
                   }
    def _record(self):
        return list(SeqRecordFactory().create(self))[0]
    record = property(_record)
    
    def __repr__(self):
        return "<Sequence id:%i len:%i '%s...'>" % (self.id, len(self.sequence), self.sequence[0:10])
    
    def biopython(self,**kwargs):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        return SeqRecord(Seq(self.sequence), str(self.id), **kwargs) 

class SequenceAc(BaseStruct, Base):
    __tablename__ = 'sequenceAc'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['sequence_key'], ['domain.domain_sequence_key']),
            {'autoload':True}
            )
    
    common_name = {'gi':'NCBI GI',
                   'ac':'Accession',
                   'ac2':'Accession2'
                   }
    
    def __repr__(self):
        return "<SequenceAc id:{0} seq:{1} prot:{2} ac:{3} ac2:{4}>".format(self.id, self.sequence_key, self.protein_key, self.ac, self.ac2) 

    def names(self):
        if self.ac and self.ac2 is None:
            return None
        return [name for name in (self.ac,self.ac2)]

class SequenceGeneName(BaseStruct, Base):
    __tablename__ = 'sequenceGeneName'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            {'autoload':True}
        )

    def __repr__(self):
        if self.db:
            return "<SequenceGeneName seq: {0} ac: {1} db: {2} symbol: {3}>".format(
                self.sequence_key, self.ac, self.db, self.symbol)
        else:
            return "<SequenceGeneName seq: {0} ac: {1} symbol: {3}>".format(
                self.sequence_key, self.ac, self.symbol)
            
class MouseIDMap(BaseStruct, Base):
    __tablename__ = 'mouseIDmap'
    __table_args__ = (
            {'autoload':True}
            )
    
    def __repr__(self):
        return "<Mouse IDs: Sequence key: {0}, GeneID: {1}, {2}, Uniprot: {3}>".format(
            self.sequence_key, self.gene_id, self.mgi_id, self.uniprot_id)
       

class MouseGeneExprCorr(BaseStruct, Base):
    __tablename__ = 'mouse_geneexpr_corr'
    __table_args__ = (
            {'autoload':True}
            )

    def __repr__(self):
        return "<Mouse gene expression data corr.: MGI Source: {0}, MGI Compared: {1}, Correlation: {2}>".format(
            self.mgi_src, self.mgi_cmp, self.correlation)

# This class represents a view. May behave differently
# class MouseFuncDomain(BaseStruct, Base):
#     __tablename__ = 'mousefunc_domains'
#     __table_args__ = (
#             {'autoload':True}
#             )
#
#     def __repr__(self):
#         return "<DomainID: {0}, ProteinID: {1}, MGI: {2}, PDB: {1}{2}>".format(self.domain_id, self.protein_id, self.mgi_id, self.pdb_id, self.chain)

# UPDATE: Based on a view. Removing.
# class Disopred(BaseStruct, Base):
#     """
#     This class represents sequecne disorder (via disopred) data, and is based on a VIEW
#     """
#     __tablename__= 'sequenceDisopred'
#     __table_args__ = (
#             ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
#             {'autoload': True},
#             )
#     id = Column(Integer, primary_key=True)
#     def __repr__(self, ):
#         return "<Disopred: sequence key {0}>".format(self.sequence_key)
# Disopred.sequence = relation(Sequence, primaryjoin="Sequence.id==Disopred.sequence_key", backref=backref("disopred",uselist=True,order_by=Disopred.ginzu_version.desc()))


class Psipred(BaseStruct, Base):
    __tablename__ = 'sequencePsiPred'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            {'autoload':True}
            )
    
    
    def _psipred(self):
        from hpf.pdb.psipred import PsipredPrediction
        return PsipredPrediction(str(self.prediction.strip()),[int(i) for i in self.confidence.strip()])
    psipred = property(_psipred)

Psipred.sequence = relation(Sequence, primaryjoin="Sequence.id==Psipred.sequence_key",
    backref=backref("psipred",uselist=True,order_by=Psipred.ginzu_version.desc()))

class PsipredFactory(object):
    def create(self,psipred,**kwargs):
        return Psipred(prediction=str(psipred),
                       confidence="".join([str(i) for i in psipred.weights]),
                       **kwargs)

class NrTaxName(BaseStruct, Base):
    __tablename__ = 'nrTaxName'
    __table_args__ = (
            {'autoload':True}
            )

class SCOPDescription(BaseStruct, Base):
    __tablename__ = 'scop_des'
    __table_args__ = (
            {'autoload':True}
            )


class SCOPClass(BaseStruct, Base):
    __tablename__ = 'scop_cla'
    __table_args__ = (
            ForeignKeyConstraint(['sccs'], ['scop_des.sccs']),
            {'autoload':True}
            )
    description = relation(SCOPDescription,
        primaryjoin="and_(SCOPDescription.sccs==SCOPClass.sccs, SCOPDescription.entrytype.in_(['cl','cf','sf','fa','dm']))",
        uselist=False)

class NrTaxNode(BaseStruct, Base):
    __tablename__ = 'nrTaxNode'
    __table_args__ = (
            ForeignKeyConstraint(['taxonomy_id'], ['nrTaxName.taxonomy_id']),
            ForeignKeyConstraint(['parent_taxonomy_id'], ['nrTaxNode.taxonomy_id']),
            ForeignKeyConstraint(['taxonomy_id'], ['nrTaxNode.parent_taxonomy_id']),
            ForeignKeyConstraint(['taxonomy_id'], ['sequenceAc.taxonomy_id']),
            {'autoload':True}
            )


    names = relation(NrTaxName, primaryjoin="NrTaxNode.taxonomy_id==NrTaxName.taxonomy_id",backref="node")
    scientific_name = relation(NrTaxName, primaryjoin="and_(NrTaxNode.taxonomy_id==NrTaxName.taxonomy_id,NrTaxName.category=='scientific name')", foreign_keys=[NrTaxName.__table__.c.taxonomy_id], uselist=False, lazy=False)
    ac = relation(SequenceAc, primaryjoin="NrTaxNode.taxonomy_id==SequenceAc.taxonomy_id",lazy=True,backref=backref("taxonomy",uselist=False))
    
NrTaxNode.parent = relation(NrTaxNode,
                            primaryjoin="NrTaxNode.taxonomy_id==NrTaxNode.parent_taxonomy_id",
                            uselist=False,
                            backref=backref("children",remote_side="NrTaxNode.parent_taxonomy_id",uselist=True))

class Experiment(BaseStruct, Base):
    __tablename__ = 'experiment'
    __table_args__ = (
            ForeignKeyConstraint(['taxonomy_id'], ['nrTaxName.taxonomy_id']),
            {'autoload':True}
            )
    
    taxon = relation(NrTaxName, primaryjoin="and_(Experiment.taxonomy_id==NrTaxName.taxonomy_id,NrTaxName.category=='scientific name')",
                     foreign_keys=[NrTaxName.__table__.c.taxonomy_id], uselist=False, lazy=False)
    
    def species(self):
        if self.taxon == None:
            return ""
        else:
            return str(self.taxon.name)
    
    def __repr__(self):
        return "<Experiment id:%i name:%s tax:%i>" % (self.id, self.name, self.taxonomy_id) 

class Protein(BaseStruct, Base):
    __tablename__ = 'protein'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['experiment_key'], ['experiment.id']),
            ForeignKeyConstraint(['sequence_key'], ['domain_sccs.parent_sequence_key']),
            ForeignKeyConstraint(['sequence_key'], ['filesystemOutfile.parent_sequence_key']),
            ForeignKeyConstraint(['sequence_key'], ['sequenceAc.sequence_key']),
            ForeignKeyConstraint(['id'], ['sequenceAc.protein_key']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    sequence = relation(Sequence, primaryjoin="Sequence.id==Protein.sequence_key", uselist=False)
    experiment = relation(Experiment, primaryjoin="Experiment.id==Protein.experiment_key",uselist=False) 
    ac = relation(SequenceAc, primaryjoin="Protein.id==SequenceAc.protein_key", uselist=False, backref=backref("protein",uselist=False))
    
    def publications(self):
        for ac in self.ac:
            for pub in ac.publications:
                yield pub
    
    def _latest_domains(self):
        # From list of all domains, get only domains from latest ginzu version.
        if not self.all_domains:
            return []

        all_domains = sorted(self.all_domains, key=lambda domain: domain.ginzu_version, reverse=True)
        latest_version = all_domains[0].ginzu_version
        latest_domains = []
        for domain in all_domains:
            if domain.ginzu_version == latest_version:
                latest_domains.append(domain)
            else:
                break
        return latest_domains
    domains = property(_latest_domains)
    
    def __repr__(self):
        return "<Protein id:%i exp:%i seq-key:%i>" % (self.id, self.experiment_key, self.sequence_key)
    
    def __json__(self):
        return {'id':self.id,
                'sequence_key':self.sequence_key,
                'experiment_key':self.experiment_key}

class Domain(BaseStruct, Base):
    __tablename__ = 'domain'
    __table_args__ = (
            ForeignKeyConstraint(['domain_sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['parent_sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['parent_sequence_key'], ['protein.sequence_key']),
            ForeignKeyConstraint(['domain_sequence_key'], ['domain_sccs.domain_sequence_key']),
            ForeignKeyConstraint(['outfile_key'], ['filesystemOutfile.id']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    sequence = relation(Sequence, primaryjoin="Sequence.id==Domain.domain_sequence_key", uselist=False)
    parent_sequence = relation(Sequence, primaryjoin="Sequence.id==Domain.parent_sequence_key", uselist=False)
    proteins = relation(Protein, primaryjoin="Protein.sequence_key==Domain.parent_sequence_key", uselist=True)

    # Cutoffs for astral finding
    ASTRAL_MCM_CUTOFF = 0.8

    def __repr__(self):
        return "<Domain id:%i seq:%i parent:%i type:%s>" % (self.id, self.domain_sequence_key, self.parent_sequence_key, self.domain_type) 

    def _known_type(self, ):
        if self.domain_type in ("psiblast", "fold_recognition"):
            return True
        else:
            return False
    known_type = property(_known_type)
    
    def _psr(self):
        if self.parent_id:
            session = object_session(self)
            return session.query(PDBSeqRes).filter_by(sequence_key=self.parent_id[3:]).all()
        else:
            return None
    pdbseqs = property(_psr)

    def _best_astral(self):
        if self.mcmdata:
            best_mcm = self.mcmdata[0]
            for mcm in self.mcmdata[1:]:
                if mcm.probability > best_mcm.probability:
                    best_mcm = mcm
            if float(best_mcm.probability) >= self.ASTRAL_MCM_CUTOFF and best_mcm.astral:
                return best_mcm.astral
        # If no conditions are met as above, return None object
        return None
    best_astral = property(_best_astral)

    # Removed, view
    # def _ginzu_version(self):
    #     if self.ginzu_run:
    #         return self.ginzu_run.ginzu_version
    #     else:
    #         return None
    # ginzu_version = property(_ginzu_version)


class DomainRegion(BaseStruct, Base):
    __tablename__ = 'domainRegion'
    __table_args__ = (
            ForeignKeyConstraint(['domain_key'], ['domain.id']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    domain = relation(Domain, primaryjoin="DomainRegion.domain_key==Domain.id", uselist=False, backref=backref("region",uselist=False))

    def __repr__(self):
        return "<DomainRegion>"# % (self.id, self.domain_sequence_key, self.parent_sequence_key, self.domain_type) 


class DomainSCCS(BaseStruct, Base):
    __tablename__ = 'domain_sccs'
    __table_args__ = (
            ForeignKeyConstraint(['domain_sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['parent_sequence_key'], ['protein.sequence_key']),
            ForeignKeyConstraint(['domain_sequence_key'], ['domain.domain_sequence_key']),
            ForeignKeyConstraint(['sccs'], ['scop_des.sccs']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    domain = relation(Domain, primaryjoin="DomainSCCS.domain_sequence_key==Domain.domain_sequence_key", uselist=False, backref=backref("sccs",uselist=False))
    protein = relation(Protein, primaryjoin="DomainSCCS.parent_sequence_key==Protein.sequence_key", uselist=False, backref=backref('sccs',uselist=True)) 
    scopDescription = relation(SCOPDescription,primaryjoin="DomainSCCS.sccs==SCOPDescription.sccs", uselist=False)

    def __repr__(self):
        return "<DomainSCCS parent_seq: %i domain_seq:%i sccs:%s conf:%f>" % (self.parent_sequence_key, self.domain_sequence_key, self.sccs, self.confidence) 

# Stores ratios on how much of astral A overlaps domain D
class AstralDomainOverlap(BaseStruct, Base):
    __tablename__ = 'astral_domain_overlap'
    __table_args__ = (
            ForeignKeyConstraint(['astral_id'], ['astral.id']),
            ForeignKeyConstraint(['astral_sid'], ['astral.sid']),
            ForeignKeyConstraint(['domain_id'], ['domain.id']),
            {'autoload':True}
            )
    domain = relation(Domain, primaryjoin="AstralDomainOverlap.domain_id==Domain.id", uselist=False, backref=backref("astral_domain_overlap", uselist=True))
    
    def __repr__(self):
        return "<AstralDomainOverlap astral: {3} ({4}, {5}-{6}), domain: {0} ({1}-{2}), pdb: {7}{8}, overlap: {9} >".format(self.domain_id, self.domain_start, self.domain_stop, self.astral_id, self.astral_sid, self.astral_start, self.astral_stop, self.pdb_id, self.chain, self.overlap)

# Stores ratios on how much of domain D overlaps astral A
class DomainAstralOverlap(BaseStruct, Base):
    __tablename__ = 'domain_astral_overlap'
    __table_args__ = (
            ForeignKeyConstraint(['domain_id'], ['domain.id']),
            ForeignKeyConstraint(['astral_id'], ['astral.id']),
            ForeignKeyConstraint(['astral_sid'], ['astral.sid']),
            {'autoload':True}
            )
    domain = relation(Domain, primaryjoin="DomainAstralOverlap.domain_id==Domain.id", uselist=False, backref=backref("domain_astral_overlap", uselist=True))
    
    def __repr__(self):
        return "<DomainAstralOverlap domain: {0} ({1}-{2}), astral: {3} ({4}, {5}-{6}), pdb: {7}{8}, overlap: {9} >".format(self.domain_id, self.domain_start, self.domain_stop, self.astral_id, self.astral_sid, self.astral_start, self.astral_stop, self.pdb_id, self.chain, self.overlap)

# Removed, from view
# class GinzuRun(BaseStruct, Base):
# # NOTE: This ORM class is based on the VIEW hpf.ginzu_run (of ddbCommon.ginzuRun)
#     __tablename__ = 'ginzu_run'
#     __table_args__ = (
#             ForeignKeyConstraint(['id'], ['domain.ginzu_key']),
#             ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
#             ForeignKeyConstraint(['sequence_key'], ['domain.domain_sequence_key']),
#             {'autoload':True}
#             )
#     id = Column(Integer, primary_key=True)
#     domain = relation(Domain, primaryjoin="GinzuRun.id==Domain.ginzu_key", uselist=False)

#     def __repr__(self):
#         return "<Ginzu run: ID {0}, sequence {1}, ginzu_version {2}>".format(self.id, self.sequence_key, self.ginzu_version)

class Structure(BaseStruct, Base):
    __tablename__ = 'structure'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['sequence_key'], ['domain.domain_sequence_key']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    
    # Relations
    domain = relation(Domain, primaryjoin="Structure.sequence_key==Domain.domain_sequence_key")
    sequence = relation(Sequence, primaryjoin="Structure.sequence_key==Sequence.id")
    
    # File storage as compressed text (string)
    compress_file_content =  deferred(Column(Binary()))
    _text = column_property(select(["uncompress(compress_file_content)"]),deferred=True)
    def _get_text(self):
        return self._text
    def _set_text(self, data):
        self.compress_file_content = func.compress(data)
        self.sha1 = func.sha1(data)
    text = synonym('_text', descriptor=property( _get_text, _set_text))
    
    def __repr__(self):
        return "<Structure id:%s seq:%s comment:%s>" % (str(self.id), str(self.sequence_key), str(self.comment))
    
    def file(self):
        """@deprecated: Use the deferred column property instead"""
        return self.text
    
    def _pdb(self):
        from Bio.PDB.PDBParser import PDBParser
        from cStringIO import StringIO
        io = StringIO()
        io.write(self.text)
        io.flush()
        io.seek(0)
        parser=PDBParser(PERMISSIVE=1)
        structure=parser.get_structure(str(self.id), io)
        return structure
    pdb = property(_pdb)

class Astral(BaseStruct, Base):
    __tablename__ = 'astral'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['structure_key'], ['structure.id']),
            ForeignKeyConstraint(['sid'], ['scop_cla.sid']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    sequence = relation(Sequence, primaryjoin="Astral.sequence_key==Sequence.id")
    structure = relation(Structure, primaryjoin="Astral.structure_key==Structure.id")
    scop = relation(SCOPClass,primaryjoin="Astral.sid==SCOPClass.sid",uselist=False)

class McmData(BaseStruct, Base):
    __tablename__ = 'mcm'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['outfile_key'], ['domain.outfile_key']),
            ForeignKeyConstraint(['structure_key'], ['structure.id']),
            ForeignKeyConstraint(['astral_structure_key'], ['astral.structure_key']),
            #ForeignKeyConstraint(['experiment_astral_ac'], ['structure.id']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    domain    = relation(Domain, primaryjoin="McmData.outfile_key==Domain.outfile_key", uselist=False, backref=backref("mcmdata",uselist=True,primaryjoin="and_(McmData.outfile_key==Domain.outfile_key,McmData.scop=='1.75')"))
    structure = relation(Structure, primaryjoin="Structure.id==McmData.structure_key", backref="mcmdata")
    astral    = relation(Astral, primaryjoin="Astral.structure_key==McmData.astral_structure_key", backref="mcmdata", uselist=False)
    protein   = relation(Protein, secondary=Domain.__table__, primaryjoin="McmData.outfile_key==Domain.outfile_key",secondaryjoin="Domain.parent_sequence_key==Protein.sequence_key", backref="mcm")
    #filesystemoutfile = relation(FilesystemOutfile, primaryjoin="McmData.sequence_key==FilesystemOutfile.sequence_key")

    def __cmp__(self, other):
        # Will sort lowest to highest. NOTE that highest MCM scores are best.
        return cmp(self.probability,other.probability)
    
    def __repr__(self):
        return "<McmData id: {0} seq: {1} sccs: {2} conf: {3}>".format(self.id, self.sequence_key, self.sccs, self.probability)

    def str_long(self, ):
        """Return a long string representation"""
        return "<McmData id: {0} seq: {5} sccs: {1} conf: {2} ratio: {3} class: {4}>".format(self.id, 
                    self.sccs, self.probability, self.ratio, self.__dict__['class'], self.sequence_key )


class Mammoth(BaseStruct, Base):
    __tablename__ = 'mammoth'
    __table_args__ = (
            ForeignKeyConstraint(['p_structure_key'], ['structure.id']),
            ForeignKeyConstraint(['e_structure_key'], ['structure.id']),
            {'autoload':True}
            )
    prediction = relation(Structure, primaryjoin="Mammoth.p_structure_key==Structure.id")
    experiment = relation(Structure, primaryjoin="Mammoth.e_structure_key==Structure.id")

    def __repr__(self):
        return "<Mammoth pred struct key: {0} exp struct key {1} zscore {2}>".format(self.p_structure_key, self.e_structure_key, self.zscore)

class FilesystemOutfile(BaseStruct, Base):
    __tablename__ = 'filesystemOutfile'
    __table_args__ = (
        ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
        ForeignKeyConstraint(['parent_sequence_key'], ['sequence.id']),
        ForeignKeyConstraint(['parent_sequence_key'], ['protein.sequence_key']),
        {'autoload':True}
    )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    sequence = relation(Sequence, primaryjoin="FilesystemOutfile.sequence_key==Sequence.id")
    protein_sequence = relation(Sequence, primaryjoin="FilesystemOutfile.parent_sequence_key==Sequence.id")
    parent_sequence = relation(Sequence, primaryjoin="FilesystemOutfile.parent_sequence_key==Sequence.id")

    def __repr__(self, ):
        return "<FilesystemOutfile: {0}, Prediction: {1}, Sequence: {2}>".format(self.id, self.prediction_code, self.sequence_key)

class DomainFoldableMap(BaseStruct, Base):
    __tablename__ = 'domain_foldable_map'
    __table_args__ = (
        ForeignKeyConstraint(['fold_sequence_key'], ['sequence.id']),
        ForeignKeyConstraint(['parent_sequence_key'], ['sequence.id']),
        ForeignKeyConstraint(['domain_sequence_key'], ['sequence.id']),
        {'autoload':True}
    )
    sequence = relation(Sequence, primaryjoin="DomainFoldableMap.fold_sequence_key==Sequence.id")
    parent_sequence = relation(Sequence, primaryjoin="DomainFoldableMap.parent_sequence_key==Sequence.id") 

class McmResultFile(BaseStruct, Base):
    __tablename__ = 'filesystemOutfileMcmResultFile'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['outfile_key'], ['filesystemOutfile.id']),
            ForeignKeyConstraint(['id'], ['mcm.outfile_mcm_result_key']),
            {'autoload':True}
            )
    sequence = relation(Sequence, primaryjoin="McmResultFile.sequence_key==Sequence.id")
    mcmdata = relation(McmData, primaryjoin="McmResultFile.id==McmData.outfile_key", backref="result_file")
    outfile = relation(FilesystemOutfile, primaryjoin="McmResultFile.outfile_key==FilesystemOutfile.id", backref="mcm_result_file")
    
    compress_file_content =  deferred(Column(Binary()))
    _text = column_property(select(["uncompress(compress_file_content)"]),deferred=True)
    def _get_text(self):
        return self._text
    def _set_text(self, data):
        self.compress_file_content = func.compress(data)
        self.sha1 = func.sha1(data)
    text = synonym('_text',
                   descriptor=property(
                                       _get_text,
                                       _set_text))
    def __repr__(self):
        return "<McmResultFile id:%i seq:%i outfile:%i>" % (self.id, self.sequence_key, self.outfile_key)

class McmRun(BaseStruct, Base):
    """Represents the hpf.mcmRun table, keeping track of foldables that have been MCMed"""
    __tablename__ = 'mcmRun'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['prediction_code'], ['domain.ibm_prediction_code']),
            ForeignKeyConstraint(['prediction_code'], ['filesystemOutfile.prediction_code']),
            {'autoload':True}
            )
    id = Column(Integer, primary_key=True)
    domain = relation(Domain, primaryjoin="McmRun.prediction_code==Domain.ibm_prediction_code", uselist=False)
    foldable = relation(FilesystemOutfile, primaryjoin="McmRun.prediction_code==FilesystemOutfile.prediction_code", uselist=False)

    def __repr__(self, ):
        return "<McmRun code: {0}, sequence: {1}, version: {2}>".format(self.prediction_code, self.sequence_key, self.version)

class RosettaCluster(BaseStruct, Base):
    __tablename__ = 'rosetta_cluster'
    __table_args__ = (
            ForeignKeyConstraint(['structure_key'], ['structure.id']),
            ForeignKeyConstraint(['convergence_key'], ['rosetta_convergence.id']),
            {'autoload':True}
            )
    structure = relation(Structure,primaryjoin="RosettaCluster.structure_key==Structure.id",lazy=True,uselist=False)

class RosettaConvergence(BaseStruct, Base):
    __tablename__ = 'rosetta_convergence'
    __table_args__ = (
            ForeignKeyConstraint(['outfile_key'], ['filesystemOutfile.id']),
            {'autoload':True}
            )

    outfile = relation(FilesystemOutfile, primaryjoin="RosettaConvergence.outfile_key==FilesystemOutfile.id", lazy=True, backref="convergence", uselist=False)
    clusters = relation(RosettaCluster, primaryjoin="RosettaConvergence.id==RosettaCluster.convergence_key", lazy=True, uselist=True,)
    cluster_centers = relation(Structure, secondary=RosettaCluster.__table__, 
                               primaryjoin="RosettaConvergence.id==RosettaCluster.convergence_key", 
                               secondaryjoin="RosettaCluster.structure_key==Structure.id", 
                               lazy=True, uselist=True,
                              )
    def __repr__(self, ):
        return "<RosettaConvergence: {0}, Prediction {1}. Radius1 {2}, Size1 {3}. Decoys {4}>".format(self.id, self.outfile.prediction_code, self.radius1, self.size1, self.total_decoys)

class FunctionPrediction(BaseStruct, Base):
    __tablename__ = 'bayes_golite_062009_3'
    __table_args__ = (
            ForeignKeyConstraint(['parent_sequence_key'], ['protein.sequence_key']),
            ForeignKeyConstraint(['domain_sequence_key'], ['domain.domain_sequence_key']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    
    def __repr__(self):
        return "<FunctionPrediction seq:%i acc:%s pls_llr:%f>" % (self.domain_sequence_key, self.mf_acc, self.pls_llr)

class Publication(BaseStruct, Base):
    __tablename__ = 'publications'
    __table_args__ = (
            ForeignKeyConstraint(['gi'], ['sequenceAc.gi']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    ac = relation(SequenceAc, primaryjoin="Publication.gi==SequenceAc.gi", uselist=False)
    
    def __repr__(self):
        return "<Publication gi:%i title:%s jour:%s date:%s>" % (self.gi, self.title[0:10], self.journal[0:10], str(self.date))

class IeaAnnotation(BaseStruct, Base):
    __tablename__ = 'hddb_iea_golite_062009'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['sequence_key'], ['protein.sequence_key']),
            ForeignKeyConstraint(['sequence_key'], ['domain.parent_sequence_key']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    proteins = relation(Protein, primaryjoin="IeaAnnotation.sequence_key==Protein.sequence_key", uselist=False, backref=backref("iea",uselist=True))
    sequence = relation(Sequence, primaryjoin="IeaAnnotation.sequence_key==Sequence.id", uselist=False)
    
    def __repr__(self):
        return "<IEA id:%i seq:%i acc:%s>" % (self.id, self.sequence_key, self.acc)

class OIDExperiment(BaseStruct, Base):
    __tablename__ = 'oid_experiments'
    __table_args__ = (
            ForeignKeyConstraint(['experiment_key'], ['experiment.id']),
            ForeignKeyConstraint(['oid_key'], ['oid.id']),
            {'autoload':True}
            )
    experiment = relation(Experiment,primaryjoin="OIDExperiment.experiment_key==Experiment.id",lazy=True)

class OID(BaseStruct, Base):
    __tablename__ = 'oid'
    __table_args__ = (
        #ForeignKeyConstraint(['id'], ['oid_experiments.oid_key']),
            {'autoload':True}
            )
    experiments = relation(Experiment,secondary=OIDExperiment.__table__,primaryjoin="OID.id==OIDExperiment.oid_key",secondaryjoin="OIDExperiment.experiment_key==Experiment.id",lazy=True,uselist=True,backref="oid")
    def _manual(self):
        return [family for family in self.families if family.manually_curated==True]
    manually_curated = property(_manual)
    def _automatic(self):
        return [family for family in self.families if family.manually_curated==False]
    oid_families = property(_automatic)


class FamilyProtein(BaseStruct, Base):
    __tablename__ = 'family_protein'
    __table_args__ = (
            ForeignKeyConstraint(['protein_key'], ['protein.id']),
            ForeignKeyConstraint(['family_key'], ['family.id']),
            {'autoload':True}
            )
    protein = relation(Protein,primaryjoin="FamilyProtein.protein_key==Protein.id",lazy=True, order_by=Protein.id)

class Family(BaseStruct, Base):
    __tablename__ = 'family'
    __table_args__ = (
            ForeignKeyConstraint(['oid_key'], ['oid.id']),
            {'autoload':True}
            )

    oid = relation(OID,primaryjoin="Family.oid_key==OID.id",lazy=True,backref=backref("families",uselist=True,lazy=True))
    proteins = relation(Protein,secondary=FamilyProtein.__table__,primaryjoin="Family.id==FamilyProtein.family_key",secondaryjoin="FamilyProtein.protein_key==Protein.id",lazy=True,uselist=True,backref="families", order_by=Protein.id)
    
    
    def __json__(self):
        return {
            'id':self.id,
            'oid_key':self.oid_key,
            'name':self.name,
            'description':self.description,
            'manually_curated':self.manually_curated
            }

    # Related familes that share proteins.
    def _related(self):
        # Aliased tables because of the double join
        fs1 = FamilyProtein.__table__.alias('fs1')
        fs2 = FamilyProtein.__table__.alias('fs2')
        return object_session(self).query(Family).join((fs1,Family.id==fs1.c.family_key)).join((fs2,fs1.c.protein_key==fs2.c.protein_key)).filter(and_(fs2.c.family_key==self.id,fs1.c.family_key!=fs2.c.family_key)).group_by(Family.id).distinct().all()
    related = property(_related)
    
    def escaped_name(self):
        return self.name.replace(" ","_")

    
    def _alignment(self):
        return self.alignments[0] if len(self.alignments)>0 else None
    alignment = property(_alignment)
    
    def _tree(self):
        return self.alignment.tree if self.alignment else None
    tree = property(_tree)
    
    def _codeml(self):
        return self.tree.codeml_models if self.tree else None
    codeml = property(_codeml)

    def _selection(self):
        if self.codeml == None:
            return []
        else:
            models = sorted(self.codeml,reverse=True,cmp= lambda x,y:cmp(x.model,y.model))
            try:
                return models[0].positive_selection
            except IndexError:
                #print "No positive selection for family {0}. Returning None".format(self.id)
                return None
    selection = property(_selection)

    def _biotree(self):
        return TreeFactory(name_func=lambda node:node.protein.ac.ac).create(self.tree.nodes, self._id)
    biotree = property(_biotree)

    def _iea(self):
        session = object_session(self)
        return session.query(IeaAnnotation).join(Family.proteins,Protein.iea).filter(Family.id==self.id).distinct().all()
    iea = property(_iea)

    def _mf(self):
        session = object_session(self)
        return session.query(IeaAnnotation).join(Family.proteins,Protein.iea).filter(Family.id==self.id).filter(IeaAnnotation.term_type=='molecular_function').distinct().all()
    mf = property(_mf)
    def _cc(self):
        session = object_session(self)
        return session.query(IeaAnnotation).join(Family.proteins,Protein.iea).filter(Family.id==self.id).filter(IeaAnnotation.term_type=='cellular_component').distinct().all()
    cc = property(_cc)
    def _bp(self):
        session = object_session(self)
        return session.query(IeaAnnotation).join(Family.proteins,Protein.iea).filter(Family.id==self.id).filter(IeaAnnotation.term_type=='biological_process').distinct().all()
    bp = property(_bp)




Family.protein_count = column_property(select([func.count(FamilyProtein.protein_key)], Family.id==FamilyProtein.family_key).correlate(Family.__table__).label('seed count'), deferred=False, group="count")

class AlignmentColcull(BaseStruct, Base):
    __tablename__ = 'alignment_colcull'
    __table_args__ = (
            ForeignKeyConstraint(['alignment_key'], ['alignment.id']),
            {'autoload':True}
            )

class AlignmentSeqcull(BaseStruct, Base):
    __tablename__ = 'alignment_seqcull'
    __table_args__ = (
            ForeignKeyConstraint(['alignment_key'], ['alignment.id']),
            {'autoload':True}
            )

class AlignmentProtein(BaseStruct, Base):
    __tablename__ = 'alignment_protein'
    __table_args__ = (
            ForeignKeyConstraint(['protein_key'], ['protein.id']),
            ForeignKeyConstraint(['alignment_key'], ['alignment.id']),

            {'autoload':True}
            )
    domains = relation(Domain,primaryjoin="AlignmentProtein.protein_key==Protein.id",secondary=Protein.__table__, secondaryjoin="Protein.sequence_key==Domain.parent_sequence_key", uselist=True,lazy=True)
    protein = relation(Protein,primaryjoin="AlignmentProtein.protein_key==Protein.id", uselist=False, lazy=True, order_by=Protein.id)
    #dna = relation(AlignmentDNA,primaryjoin="AlignmentProtein.id==AlignmentDNA.alignment_sequence_key",uselist=False,lazy=True,backref="protein")
    
class Alignment(BaseStruct, Base):
    __tablename__ = 'alignment'
    __table_args__ = (
            ForeignKeyConstraint(['family_key'], ['family.id']),
            {'autoload':True}
            )
#    compressed =  deferred(Column(Binary()))
#    text = column_property(select(["uncompress(compressed)"]),deferred=True)
    proteins = relation(AlignmentProtein,primaryjoin="Alignment.id==AlignmentProtein.alignment_key",uselist=True,lazy=True,backref=backref("alignment",uselist=False))
    culled_columns = relation(AlignmentColcull,primaryjoin="Alignment.id==AlignmentColcull.alignment_key",uselist=True,backref=backref("alignment",uselist=False))
    culled_sequences = relation(AlignmentSeqcull,primaryjoin="Alignment.id==AlignmentSeqcull.alignment_key",uselist=True,backref=backref("alignment",uselist=False))

    def _length(self):
        return len(self.proteins[0].sequence)
    length = property(_length)

    def _alignment(self):
        return AlignmentFactory().create(*self.proteins)
    alignment = property(_alignment)
Alignment.family = relation(Family,primaryjoin="Alignment.family_key==Family.id",lazy=True,backref=backref("alignments",order_by=desc(Alignment.timestamp)))
Alignment.size = column_property(select([func.count(AlignmentProtein.protein_key)], Alignment.id==AlignmentProtein.alignment_key).label('size'), deferred=False, group="count")

class TreeNode(BaseStruct, Base):
    __tablename__ = 'tree_node'
    __table_args__ = (
            ForeignKeyConstraint(['protein_key'], ['protein.id']),
            ForeignKeyConstraint(['ancestor_node'], ['tree_node.id']),
            {'autoload':True}
            )
    protein = relation(Protein,primaryjoin="TreeNode.protein_key==Protein.id",uselist=False)
TreeNode.ancestor = relation(TreeNode,remote_side=[TreeNode.__table__.c.id])

class TreeNodeDiagChars(BaseStruct, Base):
    __tablename__ = 'tree_node_diagchars'
    __table_args__ = (
            ForeignKeyConstraint(['tree_node_key'], ['tree_node.id']),
            {'autoload':True}
            )

    node = relation(TreeNode,primaryjoin="TreeNode.id==TreeNodeDiagChars.tree_node_key",uselist=False,backref=backref("diagchars",uselist=True))

class Tree(BaseStruct, Base):
    __tablename__ = 'tree'
    __table_args__ = (
            ForeignKeyConstraint(['alignment_key'], ['alignment.id']),
            ForeignKeyConstraint(['id'], ['tree_node.tree_key']),
            {'autoload':True}
            )
    
#    _descriptor = {"text":BaseStruct._compress}
#        
#    compressed =  deferred(Column(Binary()))
#    text = column_property(select(["uncompress(compressed)"]),deferred=True)
    nodes = relation(TreeNode,primaryjoin="Tree.id==TreeNode.tree_key",uselist=True,lazy=True,backref="tree",order_by=TreeNode.id.asc())
    def _tree(self):
        return TreeFactory(name_func=lambda node:node.protein.ac.ac if node.protein else None).create(self.nodes, self.id)
    tree = property(_tree)
    
    # Return the best model we have results for
    def _codeml(self):
        # Sort by Model # in reverse, ie Model 8 is the best.
        models = sorted(self.codeml_models,cmp=lambda x,y: cmp(x.model,y.model), reverse=True)
        return models[0] if len(models)>0 else None
    codeml = property(_codeml)
Tree.alignment = relation(Alignment,primaryjoin="Tree.alignment_key==Alignment.id",lazy=True,backref=backref("tree",uselist=False,order_by=Tree.timestamp.desc()))

class CodeML(BaseStruct, Base):
    __tablename__ = 'codeml'
    __table_args__ = (
            ForeignKeyConstraint(['tree_key'], ['tree.id']),
            {'autoload':True}
            )
#    compressed =  deferred(Column(Binary()))
#    text = column_property(,deferred=True)
    tree = relation(Tree, primaryjoin="CodeML.tree_key==Tree.id",lazy=True,uselist=False,backref=backref("codeml_models",uselist=True))
    
#    _descriptor = {"text":BaseStruct._compress}

class PositiveSelection(BaseStruct, Base):
    __tablename__ = 'codeml_positive_selection'
    __table_args__ = (
            ForeignKeyConstraint(['codeml_key'], ['codeml.id']),
            {'autoload':True}
            )
    def __json__(self):
        return {"column":self.column,
                "probability":float(self.probability),
                "post_mean":float(self.post_mean),
                "stderr":float(self.stderr),
                "model":self.codeml.model}
PositiveSelection.codeml = relation(CodeML, primaryjoin="PositiveSelection.codeml_key==CodeML.id", uselist=False, lazy=True,backref=backref("positive_selection",uselist=True,order_by=PositiveSelection.column.asc()))



class OntologyTerm(BaseStruct, Base):
    __tablename__ = 'mygo_term'
    __table_args__ = (
            ForeignKeyConstraint(['acc'], ['hddb_iea_golite_062009.acc']),
            ForeignKeyConstraint(['acc'], ['bayes_golite_062009_3.mf_acc']),
            {'autoload':True}
            )
    
    iea = relation(IeaAnnotation, primaryjoin="IeaAnnotation.acc==OntologyTerm.acc",backref=backref("term",uselist=False))

    def __repr__(self):
        return "<OntologyTerm acc:%s name:%s>" % (self.acc, self.name)

    def __json__(self):
        return {
            'id':self.id,
            'acc':self.acc,
            'name':self.name,
            'term_type':self.term_type
            }

class TermProbability(BaseStruct, Base):
    __tablename__ = 'probability_goLite_062009_3'
    __table_args__ = (
            ForeignKeyConstraint(['acc'], ['mygo_term.acc']),
            #ForeignKeyConstraint(['sccs'], ['scop_cla.acc']),
            ForeignKeyConstraint(['acc2'], ['mygo_term.acc']),
            {'autoload':True}
            )
    
    def __repr__(self):
        return "<TermProbability acc:%s acc2:%s probability:%f>" % (self.acc, self.acc2, self.metric)

    term1 = relation(OntologyTerm, primaryjoin="TermProbability.acc==OntologyTerm.acc",uselist=False,lazy=False)
    term2 = relation(OntologyTerm, primaryjoin="TermProbability.acc2==OntologyTerm.acc",uselist=False,lazy=False)

OntologyTerm.probability = relation(TermProbability,primaryjoin="and_(TermProbability.acc==OntologyTerm.acc, TermProbability.acc2==None)",uselist=False)

class CatalyticAtlas(BaseStruct, Base):
    __tablename__ = 'csa_2_2_11'
    __table_args__ = (
            {'autoload':True}
            )

class FireDB(Base):
    __tablename__ = 'firedb_45_full_5aug2009'
    __table_args__ = (
            {'autoload':True}
            )
    
    def __repr__(self):
        return "<FireDB pdb_chain_id:%s residues:%s>" % (self.pdb_chain_id, str(zip(self.indices(),self.resides)))

    def indices(self):
        all = []
        for i in self.residue_numbers.split(","):
            try:
                all.append(int(i))
            except:
                continue
        return all

class Matrix(BaseStruct, Base):
    __tablename__ = 'matrix'
    __table_args__ = (
            {'autoload':True}
            )
    
#    _descriptor = {"text":BaseStruct._compress}
    compressed =  deferred(Column(Binary()))
    text = BaseStruct.text #column_property(select(["uncompress(compressed)"]),deferred=True)
    def _matrix(self):
        from StringIO import StringIO
        import os
        io = StringIO()
        io.write(self.text[1:-1])
        io.flush()
        io.seek(0)
        import numpy
        return numpy.loadtxt(io)
    matrix = property(_matrix)

class Covariance(BaseStruct,Base):
    __tablename__ = 'covariance'
    __table_args__ = (
            ForeignKeyConstraint(['matrix_key'], ['matrix.id']),
            ForeignKeyConstraint(['alignment_key'], ['alignment.id']),
            ForeignKeyConstraint(['tree_key'], ['tree.id']),
            {'autoload':True}
            )
Covariance.alignment = relation(Alignment,primaryjoin=Alignment.id==Covariance.alignment_key,backref="covariance")
Covariance.zscore = relation(Matrix,primaryjoin=Matrix.id==Covariance.matrix_key)
Covariance.tree = relation(Tree,primaryjoin=Tree.id==Covariance.tree_key,backref="covariance")

class YRC(BaseStruct,Base):
    __tablename__ = 'yeastrcu'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            {'autoload':True}
            )
YRC.sequence = relation(Sequence, primaryjoin=Sequence.id==YRC.sequence_key,backref="yrc",lazy=True)

class YRCProtein(BaseStruct,Base):
    __tablename__ = 'yrc_protein'
    __table_args__ = (
            ForeignKeyConstraint(['yrc_protein_key'], ['yeastrcu.yrc_protein_key']),
            ForeignKeyConstraint(['taxonomy_id'], ['nrTaxName.taxonomy_id']),
            {'autoload':True}
            )
    yrc = relation(YRC, primaryjoin="YRC.yrc_protein_key==YRCProtein.yrc_protein_key",backref=backref("proteins",uselist=True),lazy=True)
    taxon = relation(NrTaxName, primaryjoin="and_(YRCProtein.taxonomy_id==NrTaxName.taxonomy_id,NrTaxName.category=='scientific name')",
                     foreign_keys=[NrTaxName.__table__.c.taxonomy_id], 
                     uselist=False, lazy=False)

    def __json__(self):
        return {'id':self.id,
                'yrc_protein_key':self.yrc_protein_key,
                'taxonomy_id':self.taxonomy_id,
                'name':self.name,
                'description':self.description}


class PDBIndex(BaseStruct, Base):
    __tablename__ = 'pdbIndex'
    __table_args__ = (
            ForeignKeyConstraint(['pdbId'], ['domain_sccs.pdbid']),
            {'autoload':True}
            )

DomainSCCS.pdb = relation(PDBIndex,primaryjoin="DomainSCCS.pdbid==PDBIndex.pdbId",uselist=False)

class PDBSeqRes(BaseStruct, Base):
    __tablename__ = 'pdbSeqRes'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['structure_key'], ['structure.id']),
            ForeignKeyConstraint(['pdb_key'], ['pdbIndex.id']),
            {'autoload':True}
            )
    structure = relation(Structure, primaryjoin="PDBSeqRes.structure_key==Structure.id")
    sequence = relation(Sequence, primaryjoin="PDBSeqRes.sequence_key==Sequence.id")
    pdb = relation(PDBIndex, primaryjoin="PDBSeqRes.pdb_key==PDBIndex.id",backref="seqres")
    
    def __repr__(self):
        return "<PDBSeqRes id: {0}, pdb_id: {1}, chain: {2}, sequence_key: {3}>".format(self.id, self.pdb.pdbId, self.chain, self.sequence_key) 

class StructureMammoth(BaseStruct, Base):
    __tablename__ = 'structure_mammoth'
    __table_args__ = (
            ForeignKeyConstraint(['prediction_id'], ['structure.id']),
            ForeignKeyConstraint(['experiment_id'], ['structure.id']),
            {'autoload':True}
            )

    prediction_structure = relation(Structure, primaryjoin="StructureMammoth.prediction_id==Structure.id", uselist=False)
    experiment_structure = relation(Structure, primaryjoin="StructureMammoth.experiment_id==Structure.id", uselist=False)

    def __repr__(self, ):
        return "<StructureMammoth: pred {0}, exp {1}, similarity {2}>".format(self.prediction_id, self.experiment_id, self.zscore)

class StructureMammoth01(BaseStruct, Base):
    __tablename__ = 'structure_mammoth_01'
    __table_args__ = (
            ForeignKeyConstraint(['prediction_id'], ['structure.id']),
            ForeignKeyConstraint(['experiment_id'], ['structure.id']),
            {'autoload':True}
            )

    prediction_structure = relation(Structure, primaryjoin="StructureMammoth01.prediction_id==Structure.id", uselist=False)
    experiment_structure = relation(Structure, primaryjoin="StructureMammoth01.experiment_id==Structure.id", uselist=False)

    def __repr__(self, ):
        return "<StructureMammoth01: pred {0}, exp {1}, similarity {2}>".format(self.prediction_id, self.experiment_id, self.zscore)


class StructureMammoth1186(BaseStruct, Base):
    __tablename__ = 'structure_mammoth_1186'
    __table_args__ = (
            ForeignKeyConstraint(['prediction_id'], ['structure.id']),
            ForeignKeyConstraint(['experiment_id'], ['structure.id']),
            {'autoload':True}
            )

    prediction_structure = relation(Structure, primaryjoin="StructureMammoth1186.prediction_id==Structure.id", uselist=False)
    experiment_structure = relation(Structure, primaryjoin="StructureMammoth1186.experiment_id==Structure.id", uselist=False)

    def __repr__(self, ):
        return "<StructureMammoth1186: pred {0}, exp {1}, similarity {2}>".format(self.prediction_id, self.experiment_id, self.zscore)


#kdrew: add to lookup dictionary the table names of structure_mammoth tables
#dpb: this aggravates me.
StructureMammoth_class_lookup = { 
    'default':StructureMammoth,
    'part01':StructureMammoth01,
    '1186':StructureMammoth1186, 
}

class MammothRun(BaseStruct, Base):
    """Represents the hpf.mammothRun table, keeping track of structure supergroups that have been Mammoth-ed"""
    __tablename__ = 'mammothRun'
    __table_args__ = (
            ForeignKeyConstraint(['supergroup_key'], ['mammothSupergroup.id']),
            {'autoload':True}
            )
    id = Column(Integer, primary_key=True)

    def __repr__(self, ):
        return "<MammothRun supergroup_key: {0}, group_key1: {1}, group_key2: {2}, version: {3}>".format(self.supergroup_key, self.group_key1, self.group_key2, self.version)

class MammothGroup(BaseStruct, Base):
    """Represents the hpf.mammothGroup table, keeping track of structure_keys, their supergroups and groups """
    __tablename__ = 'mammothGroup'
    __table_args__ = (
            ForeignKeyConstraint(['structure_key'], ['structure.id']),
            ForeignKeyConstraint(['supergroup_key'], ['mammothSupergroup.id']),
            {'autoload':True}
            )
    id = Column(Integer, primary_key=True)

    def __repr__(self, ):
        return "<MammothGroup supergroup_key: {0}, group_key: {1}, structure_key: {2}>".format(self.supergroup_key, self.group_key, self.structure_key)

class MammothSupergroup(BaseStruct, Base):
    """Represents the hpf.mammothSupergroup table, keeping track of supergroups """
    __tablename__ = 'mammothSupergroup'
    __table_args__ = (
            ForeignKeyConstraint(['id'], ['mammothGroup.supergroup_key']),
            ForeignKeyConstraint(['id'], ['mammothRun.supergroup_key']),
            {'autoload':True}
            )
    id = Column(Integer, primary_key=True)

    def __repr__(self, ):
        return "<MammothSupergroup id: {0}, name: {1}, comment: {2}>".format(self.id, self.name, self.comment)




class AstralComparison(BaseStruct, Base):
    __tablename__ = 'astral_mammoth'
    __table_args__ = (
            ForeignKeyConstraint(['prediction'], ['astral.sid']),
            ForeignKeyConstraint(['experiment'], ['astral.sid']),
            {'autoload':True}
            )

    source_astral = relation(Astral, primaryjoin="AstralComparison.prediction==Astral.sid", uselist=False)
    comparison_astral = relation(Astral, primaryjoin="AstralComparison.experiment==Astral.sid", uselist=False)
    
    def __init__(self, ):
        self.similarity = self.zscore

    def __repr__(self, ):
        return "<AstralComparison source: {0}, target: {1}, similarity: {2}>".format(self.prediction, self.experiment, self.zscore)

class HomologStructAllVAll(BaseStruct, Base):
    __tablename__ = 'homolog_struct_allvall_best100'
    __table_args__ = ({'autoload':True})

    def __repr__(self, ):
        return "<HomologStructAllVAll: source {0}, target {1}, score {2}>".format(self.source_id, self.target_id, self.score)


## HHPred factories

class HHPredRFFactory(object):
    """Creates an HHPredResultFile ORM object from the fields of a given hpf.hhpred.HHPredResultFile object
    A sequence and version key must be provided OR must exist as an instance variables in the hhpred_resultfile object
    Should use an hpf.hhpred.HpfHHPredResultFile object (with sequence_key and version_key).
    """
    def create(self, hhpRF, version_key=None, sequence_key=None, debug=False):
        # Check for version and sequence keys (required for DB linkery)
        if not version_key:
            try:
                version_key = hhpRF.version_key
            except AttributeError:
                print "Must provide version key of corresponding hhpred run for ORM object instantiation"
                raise
        if not sequence_key:
            try:
                sequence_key = hhpRF.sequence_key
            except AttributeError:
                print "Must provide sequence key for ORM object instantiation"
                raise
        
        # Get results file text to store (compressed) in DB
        with open(hhpRF.results_file) as handle:
            result_text = handle.read()
        if result_text == None or result_text == '':
            raise Exception("Results file '{0}' from HHpredResultFile object is empty".format(hhpRF.results_file))
        
        # Create and return ORM object
        hhprf_dbo = HHPredResultFile(sequence_key=sequence_key, 
                                version_key=version_key, 
                                filename=hhpRF.results_file, 
                                text=result_text, 
                                insert_date=datetime.now()
                                )
        if debug: print "Created {0}".format(hhprf_dbo)
        return hhprf_dbo

class HHPredHitFactory(object):
    """Creates an HHPredHit ORM object from the fields of a given hpf.hhpred.HHPredHit object.
    Must pass in resultfile key to properly link the HHPredHit in the database. Note that this
    key should be taken from the HHPredResultFile ORM object, as HHPredResultFile.id, which will
    only be set if the HHPredResultFile object has been pushed to the DB.
    """
    def create(self, hhpHit, resultfile_key=None, debug=False):
        if not resultfile_key:
            raise Exception("HHPredHit must be associated with a parent HHPredResultFile. Must provide a resultfile key.")
        hhphit_dbo = HHPredHit(resultfile_key=resultfile_key, 
                               number=hhpHit.number,
                               hit_id=hhpHit.hit_id,
                               hit_desc=hhpHit.hit_description,
                               probability=hhpHit.probability,
                               evalue=hhpHit.evalue,
                               pvalue=hhpHit.pvalue,
                               score=hhpHit.score,
                               ss_score=hhpHit.ss,
                               columns=hhpHit.columns,
                               match_states=hhpHit.match_states,
                               query_start=hhpHit.query_start,
                               query_stop=hhpHit.query_stop,
                               template_start=hhpHit.template_start,
                               template_stop=hhpHit.template_stop,
                               )
        if debug: print "Created {0}".format(hhphit_dbo)
        return hhphit_dbo

## HHPRED classes

class HHPredVersion(BaseStruct, Base):
    __tablename__ = 'hhpred_version'
    __table_args__ = ({'autoload':True})

    def __repr__(self, ):
        return "<HHPredVersion {0}: hhpred {1} run against {2} ({3}...)>".format(self.id, self.version, self.comparison_db, self.comment[:25])

class HHPredResultFile(BaseStruct, Base):
    __tablename__ = 'hhpred_result_file'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['sequence_key'], ['domain.parent_sequence_key']),
            ForeignKeyConstraint(['version_key'], ['hhpred_version.id']),
            {'autoload':True}
            )
    # Relations
    sequence = relation(Sequence, primaryjoin="HHPredResultFile.sequence_key==Sequence.id")
    version  = relation(HHPredVersion, primaryjoin="HHPredResultFile.version_key==HHPredVersion.id", backref=backref("result_files", uselist=True))

    # File storage as compressed text (string)
    compress_file_content =  deferred(Column(Binary()))
    _text = column_property(select(["uncompress(compress_file_content)"]),deferred=True)
    def _get_text(self):
        return self._text
    def _set_text(self, data):
        self.compress_file_content = func.compress(data)
        self.sha1 = func.sha1(data)
    text = synonym('_text', descriptor=property( _get_text, _set_text))

    def __repr__(self, ):
        return "<HHpredResultFile id: {3} version: {0} sequence: {1} file: {2}>".format(self.version_key, self.sequence_key, self.filename, self.id)

class HHPredHit(BaseStruct, Base):
    __tablename__ = 'sequenceHHPred'
    __table_args__ = (
            ForeignKeyConstraint(['resultfile_key'], ['hhpred_result_file.id']),
            {'autoload':True}
            )
    # Relations
    result_file = relation(HHPredResultFile, primaryjoin="HHPredHit.resultfile_key==HHPredResultFile.id", backref=backref("hits", uselist=True))
    
    def __repr__(self, ):
        return "<HHPredHit id: {0}, hit: {1} eval: {2} from resultfile: {3}>".format(self.id, self.hit_id, self.evalue, self.resultfile_key)


# Mappings that need to be made after classes are prepared
Sequence.ac = relation(SequenceAc, primaryjoin=Sequence.id==SequenceAc.sequence_key, uselist=True)
Sequence.proteins = relation(Protein, primaryjoin=Sequence.id==Protein.sequence_key, uselist=True)
Sequence.domains = relation(Domain, primaryjoin=Sequence.id==Domain.domain_sequence_key, uselist=True)
Sequence.iea = relation(IeaAnnotation, primaryjoin=Sequence.id==IeaAnnotation.sequence_key)

SequenceAc.publications = relation(Publication,primaryjoin=SequenceAc.gi==Publication.gi, uselist=True)

Experiment.protein_count = column_property(select([func.count(Protein.id)], Experiment.id==Protein.experiment_key).label('protein_count'), deferred=True, group="count")
Experiment.domain_count = column_property(select([func.count(Domain.id)], and_(Experiment.id==Protein.experiment_key,Domain.parent_sequence_key==Protein.sequence_key)).label('domain_count'), deferred=True, group="count")

Protein.all_domains = relation(Domain, primaryjoin=Protein.sequence_key==Domain.parent_sequence_key, uselist=True)
Protein.function_predictions = relation(FunctionPrediction, primaryjoin=Protein.sequence_key==FunctionPrediction.parent_sequence_key, uselist=True)
Protein.fold_recognition = column_property(select([func.count(DomainSCCS.size)>0], and_(Protein.sequence_key==DomainSCCS.parent_sequence_key, DomainSCCS.domain_type=='fold_recognition')).correlate(Protein.__table__),group="sccs", deferred=False)
Protein.rosetta = column_property(select([func.max(DomainSCCS.confidence)], and_(Protein.sequence_key==DomainSCCS.parent_sequence_key, and_(DomainSCCS.domain_type!='psiblast',DomainSCCS.domain_type!='fold_recognition'))).correlate(Protein.__table__),group="sccs", deferred=False)
Protein.psiblast = column_property(select([func.count(DomainSCCS.size)>0], and_(Protein.sequence_key==DomainSCCS.parent_sequence_key, DomainSCCS.domain_type=='psiblast')).correlate(Protein.__table__), group="sccs", deferred=False)
Protein.experiment_name = column_property(select([Experiment.name],Protein.experiment_key==Experiment.id).correlate(Protein.__table__), group="tax", deferred=False)
Protein.organism = column_property(select([NrTaxName.name],and_(Protein.experiment_key==Experiment.id, Experiment.taxonomy_id==NrTaxName.taxonomy_id, NrTaxName.category=='scientific name')).correlate(Protein.__table__),group="tax", deferred=False)

Domain.function_predictions = relation(FunctionPrediction, primaryjoin=Domain.domain_sequence_key==FunctionPrediction.domain_sequence_key, uselist=True)
Domain.iea = relation(IeaAnnotation, primaryjoin=Domain.parent_sequence_key==IeaAnnotation.sequence_key)
Domain.experiments = relation(Experiment, secondary=Protein.__table__, primaryjoin="Domain.parent_sequence_key==Protein.sequence_key",secondaryjoin="Protein.experiment_key==Experiment.id", uselist=True)
Domain.structures = relation(Structure, primaryjoin="Structure.sequence_key==Domain.domain_sequence_key", uselist=True)
# Removed, from view
# Domain.ginzu_run = relation(GinzuRun, primaryjoin="Domain.ginzu_key==GinzuRun.id", uselist=False)
Domain.foldable = relation(FilesystemOutfile, primaryjoin="FilesystemOutfile.id==Domain.outfile_key", uselist=False)

DomainAstralOverlap.astral = relation(Astral, primaryjoin="DomainAstralOverlap.astral_id==Astral.id", uselist=False, backref=backref("domain_astral_overlap", uselist=True))
AstralDomainOverlap.astral = relation(Astral, primaryjoin="AstralDomainOverlap.astral_id==Astral.id", uselist=False, backref=backref("astral_domain_overlap", uselist=True))

FunctionPrediction.p_term = relation(OntologyTerm, primaryjoin="FunctionPrediction.p_acc==OntologyTerm.acc",uselist=False,lazy=False)
FunctionPrediction.l_term = relation(OntologyTerm, primaryjoin="FunctionPrediction.l_acc==OntologyTerm.acc",uselist=False,lazy=False)


#Family.sequences = relation(Sequence,secondary=FamilySequence.__table__,primaryjoin="Family.id==FamilySequence.family_key",secondaryjoin="FamilySequence.sequence_key==Sequence.id")
#Family.seed_proteins = relation(Protein,secondary=FamilySequence.__table__,primaryjoin="Family.id==FamilySequence.family_key",secondaryjoin="and_(FamilySequence.sequence_key==Protein.sequence_key, Protein.experiment_key==1162)")
#Family.hpf_proteins = relation(Protein,secondary=FamilySequence.__table__,primaryjoin="Family.id==FamilySequence.family_key",secondaryjoin="and_(FamilySequence.sequence_key==Protein.sequence_key, Protein.experiment_key!=1162)")
#Family.sequence_count = column_property(select([func.count(Sequence.id)], and_(Family.id==FamilySequence.family_key,FamilySequence.sequence_key==Sequence.id)).label('sequence_count'))
#Family.seed_count = column_property(select([func.count(Protein.id)], and_(Family.id==FamilySequence.family_key,FamilySequence.sequence_key==Protein.sequence_key,Protein.experiment_key==1162)).label('seed_count'))
#Family.hpf_count = column_property(select([func.count(Protein.id)], and_(Family.id==FamilySequence.family_key,FamilySequence.sequence_key==Protein.sequence_key,Protein.experiment_key!=1162)).label('hpf_count'))
#Family.tree = relation(Tree,primaryjoin=Family.id==Tree.family_key,lazy=True, uselist=False, order_by=Tree.timestamp.desc())
#Family.alignment = relation(FamilySequence, primaryjoin=and_(Family.id==FamilySequence.family_key, FamilySequence.seed==True, FamilySequence.alignment!=None), lazy=True, uselist=True)


OntologyTerm.probability = relation(TermProbability,primaryjoin=and_(OntologyTerm.acc==TermProbability.acc, TermProbability.acc2==None),uselist=False,lazy=False)
OntologyTerm.superfamilies = relation(TermProbability,primaryjoin=and_(OntologyTerm.acc==TermProbability.acc2,TermProbability.acc.like('%.%.%')),lazy=True,uselist=True,order_by=desc(TermProbability.metric))
OntologyTerm.common = relation(TermProbability,primaryjoin=and_(OntologyTerm.acc==TermProbability.acc2,TermProbability.acc.like('GO:%')),lazy=True,uselist=True,order_by=desc(TermProbability.metric))
                                                   #addresses_table.c.user_id==users_table.c.user_id

#McmData.protein_fs = relation(Protein, secondary=FilesystemOutfile.__table__, primaryjoin="McmData.sequence_key==FilesystemOutfile.sequence_key",secondaryjoin="FilesystemOutfile.parent_sequence_key==Protein.sequence_key", backref="mcm")

#binds = {#Domain: ddbCommon,
#         #Domain_Sccs: hpf,
#         #Experiment: hpf,
#         Protein: bddb,
#         Sequence: ddbCommon
#        }

#binds=engine,twophase=True)

        
    
