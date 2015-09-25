'''
Load and map firedb catalytic residues
Created on Jan 27, 2010
@author: Patrick
'''
'''
Created on Jan 5, 2010

@author: patrick
'''
import sqlalchemy.ext.declarative
from sqlalchemy import create_engine, ForeignKeyConstraint
from sqlalchemy.orm import sessionmaker, relation, dynamic_loader, backref, column_property
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.sql import and_, or_, select, func
from hpf.runtime import runtime
from sqlalchemy.sql.expression import desc

url = "mysql://patrick:patrick_nyu@127.0.0.1:3306/catalytic_site_atlas"
    
def clear():
    global engine, Session, Base
    runtime().debug("Clearing engine/session mapping")
    engine, Session, Base = (None,None,None)

def setup():
    e = create_engine(url)
    b = sqlalchemy.ext.declarative.declarative_base(bind=e)
    s = sessionmaker(bind=e)
    return e,b,s

engine, Base, Session = setup()

class CatalyticResidue(Base):
    __tablename__ = 'CSA_2_2_11'
    __table_args__ = (
            {'autoload':True}
            )
    
    def __repr__(self):
        return "<CatalyticResidue pdb_id:%s chain:%s residues:%i>" % (self.pdb_id, self.chain, self.residue_number)

        