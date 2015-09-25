'''
Essentially a cut-down version of Patrick Winters' db.py, used for testing limited functionality.
Feel free to edit and break in any way desired, but only do so in dev and test databases
(eg mysql:handbanana:hpf_dev)

@author Duncan Penfold-Brown, 2/17/2011
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

#url = "mysql://dpb:dpb_nyu@ms2.bio.nyu.edu:3306"
url = "mysql://dpb:dpb_nyu@handbanana.bio.nyu.edu:3306/"
#url = "mysql://kdrew:bonneau@127.0.0.1:13307/"
#url = "mysql://pfp:pfp_nyu@db/"

temp = None
engine = None
Base = None
Session = None
def clear():
    global engine, Session, Base
    runtime().debug("Clearing HDDB engine/session mapping")
    engine, Session, Base = (None,None,None)

def setup(db='hpf_dev'):
    # create_engine params: implicit_returning=True, pool_recycle=300, echo=True | echo="debug", pool_timeout=10)
    e = create_engine(url+db, echo=False) 
    b = sqlalchemy.ext.declarative.declarative_base(bind=e)
    s = sessionmaker(bind=e)
    return e,b,s

def rebind(**kwargs):
    global engine, Base, Session
    clear()
    engine, Base, Session = setup(**kwargs)

from hpf.hddb import tunnel
engine, Base, Session = setup()

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


class Sequence(BaseStruct, Base):
    __tablename__ = 'sequence'
    __table_args__ = (
            ForeignKeyConstraint(['id'], ['sequenceAc.sequence_key']),
            {'autoload':True}
            )
    
    common_name = {'sequence':'Protein Sequence'
                   }
    
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
            {'autoload':True}
            )
    
    common_name = {'gi':'NCBI GI',
                   'ac':'Accession',
                   'ac2':'Accession2'
                   }
    
    def __repr__(self):
        return "<SequenceAc id:%i seq:%i gi:%i ac:%s ac2:%s>" % (self.id, self.sequence_key, self.gi, self.ac, self.ac2) 

    def names(self):
        if self.ac and self.ac2 is None:
            return None
        return [name for name in (self.ac,self.ac2)]

class Experiment(BaseStruct, Base):
    __tablename__ = 'experiment'
    __table_args__ = (
            {'autoload':True}
            )
    
    def __repr__(self):
        return "<Experiment id:%i name:%s tax:%i>" % (self.id, self.name, self.taxonomy_id) 


class Protein(BaseStruct, Base):
    __tablename__ = 'protein'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['experiment_key'], ['experiment.id']),
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
    
    def __repr__(self):
        if self.id:
            return "<Protein id:%i exp:%i key:%i>" % (self.id, self.experiment_key, self.sequence_key)
        else:
            return "<Protein id:(Not yet reconciled) exp:%i key:%i>" % (self.experiment_key, self.sequence_key)
    
    def __json__(self):
        return {'id':self.id,
                'sequence_key':self.sequence_key,
                'experiment_key':self.experiment_key}

class MouseIDMap(BaseStruct, Base):
    __tablename__ = 'mouseIDmap'
    __table_args__ = (
            {'autoload':True}
            )
    
    def __repr__(self):
        return "<Mouse IDs: Sequence key: {0}, GeneID: {1}, {2}, Uniprot: {3}>".format(self.sequence_key, self.gene_id, self.mgi_id, self.uniprot_id)
        
        
##
## Mappings that need to be made after classes are prepared
##

Sequence.ac = relation(SequenceAc, primaryjoin=Sequence.id==SequenceAc.sequence_key, uselist=True)
Sequence.proteins = relation(Protein, primaryjoin=Sequence.id==Protein.sequence_key, uselist=True)

Experiment.protein_count = column_property(select([func.count(Protein.id)], Experiment.id==Protein.experiment_key).label('protein_count'), deferred=True, group="count")

Protein.experiment_name = column_property(select([Experiment.name],Protein.experiment_key==Experiment.id).correlate(Protein.__table__), group="tax", deferred=False)




