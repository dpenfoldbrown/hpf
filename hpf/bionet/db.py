'''
Created on Jan 11, 2010

@author: patrick
'''
assert False, "This is definitely not working yet"
import sqlalchemy.ext.declarative
from sqlalchemy import create_engine, ForeignKeyConstraint
from sqlalchemy.orm import sessionmaker, relation

url = "mysql://bionet:bionet_nyu@err.bio.nyu.edu/"
attributes_engine = create_engine(url+'hpf')
Base = sqlalchemy.ext.declarative.declarative_base(bind=engine)

class Sequence(Base):
    __tablename__ = 'sequence'
    __table_args__ = (
            ForeignKeyConstraint(['id'], ['protein.sequence_key']),
            ForeignKeyConstraint(['id'], ['domain.domain_sequence_key']),
            ForeignKeyConstraint(['id'], ['sequenceAc.sequence_key']),
            ForeignKeyConstraint(['id'], ['hddb_iea_golite_062009.sequence_key']),
            {'autoload':True}
            )
    def __repr__(self):
        return "<Sequence id:%i len:%i '%s...'>" % (self.id, len(self.sequence), self.sequence[0:10]) 

class SequenceAc(Base):
    __tablename__ = 'sequenceAc'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['sequence_key'], ['protein.sequence_key']),
            ForeignKeyConstraint(['sequence_key'], ['domain.domain_sequence_key']),
            {'autoload':True}
            )
    def __repr__(self):
        return "<SequenceAc id:%i seq:%i gi:%i ac:%s ac2:%s>" % (self.id, self.sequence_key, self.gi, self.ac, self.ac2) 

    def names(self):
        if self.ac and self.ac2 is None:
            return None
        return [name for name in (self.ac,self.ac2)]

class Experiment(Base):
    __tablename__ = 'experiment'
    __table_args__ = (
            ForeignKeyConstraint(['id'], ['protein.experiment_key']),
            {'autoload':True}
            )
    def __repr__(self):
        return "<Experiment id:%i name:%s tax:%i>" % (self.id, self.name, self.taxonomy_id) 

class Protein(Base):
    __tablename__ = 'protein'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['experiment_key'], ['experiment.id']),
            ForeignKeyConstraint(['sequence_key'], ['domain.parent_sequence_key']),
            ForeignKeyConstraint(['sequence_key'], ['domain_sccs.parent_sequence_key']),
            ForeignKeyConstraint(['sequence_key'], ['sequenceAc.sequence_key']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    sequence = relation(Sequence, primaryjoin="Sequence.id==Protein.sequence_key", uselist=False)
    experiment = relation(Experiment, primaryjoin="Experiment.id==Protein.experiment_key",uselist=False) 
    ac = relation(SequenceAc, primaryjoin="Protein.sequence_key==SequenceAc.sequence_key", uselist=True) 
    
    def publications(self):
        for ac in self.ac:
            for pub in ac.publications:
                yield pub
    
    def __repr__(self):
        return "<Protein id:%i exp:%i key:%i>" % (self.id, self.experiment_key, self.sequence_key) 
 
class Domain(Base):
    __tablename__ = 'domain'
    __table_args__ = (
            ForeignKeyConstraint(['domain_sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['parent_sequence_key'], ['protein.sequence_key']),
            ForeignKeyConstraint(['domain_sequence_key'], ['domain_sccs.domain_sequence_key']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    sequence = relation(Sequence, primaryjoin="Sequence.id==Domain.domain_sequence_key", uselist=False)
    proteins = relation(Protein, primaryjoin="Protein.sequence_key==Domain.parent_sequence_key", uselist=True)

    def __repr__(self):
        return "<Domain id:%i seq:%i parent:%i type:%s>" % (self.id, self.domain_sequence_key, self.parent_sequence_key, self.domain_type) 

class DomainSCCS(Base):
    __tablename__ = 'domain_sccs'
    __table_args__ = (
            ForeignKeyConstraint(['domain_sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['parent_sequence_key'], ['protein.sequence_key']),
            ForeignKeyConstraint(['domain_sequence_key'], ['domain.domain_sequence_key']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    domain = relation(Domain, primaryjoin="DomainSCCS.domain_sequence_key==Domain.domain_sequence_key")
    protein = relation(Protein, primaryjoin="DomainSCCS.parent_sequence_key==Protein.sequence_key") 

    def __repr__(self):
        return "<DomainSCCS id:%i seq:%i parent:%i sccs:%s conf:%f>" % (self.id, self.domain_sequence_key, self.parent_sequence_key, self.sccs, self.confidence) 

class McmData(Base):
    __tablename__ = 'mcmData'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['sequence_key'], ['domain.domain_sequence_key']),
            ForeignKeyConstraint(['structure_key'], ['structure.id']),
            #ForeignKeyConstraint(['experiment_astral_ac'], ['structure.id']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    domain = relation(Domain, primaryjoin="McmData.sequence_key==Domain.domain_sequence_key")
    
    def __repr__(self):
        return "<McmData id:%i seq:%i sccs:%s conf:%f>" % (self.id, self.sequence_key, self.experiment_sccs, self.probability) 

class Structure(Base):
    __tablename__ = 'structure'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['sequence_key'], ['domain.domain_sequence_key']),
            ForeignKeyConstraint(['id'], ['domain.structure_key']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    domain = relation(Domain, primaryjoin="Structure.sequence_key==Domain.domain_sequence_key")
    sequence = relation(Sequence, primaryjoin="Structure.sequence_key==Sequence.id")
    def __repr__(self):
        return "<Structure id:%i seq:%i>" % (self.id, self.sequence_key)
    
    def file(self):
        return engine.execute("select uncompress(compress_file_content) as file from structure where id=%i" % self.id).fetchone()

class FunctionPrediction(Base):
    __tablename__ = 'bayes_golite_062009_3'
    __table_args__ = (
            ForeignKeyConstraint(['parent_sequence_key'], ['protein.sequence_key']),
            ForeignKeyConstraint(['domain_sequence_key'], ['domain.domain_sequence_key']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    
    def __repr__(self):
        return "<FunctionPrediction seq:%i acc:%s pls_llr:%f>" % (self.domain_sequence_key, self.mf_acc, self.pls_llr)

class Publication(Base):
    __tablename__ = 'publications'
    __table_args__ = (
            ForeignKeyConstraint(['gi'], ['sequenceAc.gi']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    ac = relation(SequenceAc, primaryjoin="Publication.gi==SequenceAc.gi", uselist=False)
    
    def __repr__(self):
        return "<Publication gi:%i title:%s jour:%s date:%s>" % (self.gi, self.title[0:10], self.journal[0:10], str(self.date))

class IeaAnnotation(Base):
    __tablename__ = 'hddb_iea_golite_062009'
    __table_args__ = (
            ForeignKeyConstraint(['sequence_key'], ['sequence.id']),
            ForeignKeyConstraint(['sequence_key'], ['protein.sequence_key']),
            ForeignKeyConstraint(['sequence_key'], ['domain.parent_sequence_key']),
            {'autoload':True}
            )
    #__mapper_args__ = {"properties":{"sequence":relation(Sequence,primaryjoin="Sequence.id==Protein.sequence_key")}}
    protein = relation(Protein, primaryjoin="IeaAnnotation.sequence_key==Protein.sequence_key", uselist=False)
    sequence = relation(Sequence, primaryjoin="IeaAnnotation.sequence_key==Sequence.id", uselist=False)
    
    def __repr__(self):
        return "<IEA id:%i seq:%i acc:%s>" % (self.id, self.sequence_key, self.acc)

# Mappings that need to be made after classes are prepared
Sequence.ac = relation(SequenceAc, primaryjoin=Sequence.id==SequenceAc.sequence_key, uselist=True)
Sequence.proteins = relation(Protein, primaryjoin=Sequence.id==Protein.sequence_key, uselist=True)
Sequence.domains = relation(Domain, primaryjoin=Sequence.id==Domain.domain_sequence_key, uselist=True)
Sequence.iea = relation(IeaAnnotation, primaryjoin=Sequence.id==IeaAnnotation.sequence_key)

SequenceAc.publications = relation(Publication,primaryjoin=SequenceAc.gi==Publication.gi, uselist=True)

Experiment.proteins = relation(Protein, primaryjoin=Experiment.id==Protein.experiment_key, uselist=True)

Protein.domains = relation(Domain, primaryjoin=Protein.sequence_key==Domain.parent_sequence_key, uselist=True)
Protein.sccs = relation(DomainSCCS, primaryjoin=Protein.sequence_key==DomainSCCS.parent_sequence_key, uselist=True)
Protein.function_predictions = relation(FunctionPrediction, primaryjoin=Protein.sequence_key==FunctionPrediction.parent_sequence_key, uselist=True)
Protein.iea = relation(IeaAnnotation, primaryjoin=Protein.sequence_key==IeaAnnotation.sequence_key)

Domain.sccs = relation(DomainSCCS, primaryjoin=DomainSCCS.domain_sequence_key==Domain.domain_sequence_key, uselist=False)
Domain.function_predictions = relation(FunctionPrediction, primaryjoin=Domain.domain_sequence_key==FunctionPrediction.domain_sequence_key, uselist=True)
Domain.iea = relation(IeaAnnotation, primaryjoin=Domain.parent_sequence_key==IeaAnnotation.sequence_key)
Domain.mcmdata = relation(McmData, primaryjoin=Domain.domain_sequence_key==McmData.sequence_key)

McmData.structure = relation(Structure, primaryjoin=McmData.structure_key==Structure.id)

#binds = {#Domain: ddbCommon,
#         #Domain_Sccs: hpf,
#         #Experiment: hpf,
#         Protein: bddb,
#         Sequence: ddbCommon
#        }

Session = sessionmaker()#binds=engine,twophase=True)
