import sys
from hpf.xml import SAXFactory, Attributes
from xml.sax import saxutils
from collections import defaultdict

from hpf.mapper import ProteinDomainMapper, ProteinAlignmentMapper
from hpf.pdb import get_crystal
from hpf.pdb.mapper import ProteinToPDBAtom
from hpf.pdb.dssp import SHORT_NAME, SS, ASA, RASA, DSSPFeatureSet

class Structure(object):
    """
    Simple struct for the indexing mechanism for structure types.
    """
    def __init__(self, id, chain=None):
        self.index = str(id)+str(chain) if chain else str(id)
        self.id = id
        self.chain = chain
        self.type = "pdb" if chain else "decoy"

    def __eq__(self, other):
        return self.index == other.index

    def __hash__(self):
        return hash(self.index)

    def __cmp__(self, other):
        return cmp(self.index, other.index)

class ColumnStructureAttribute(object):
    """
    Encapsulates all of the required information for a column's mapping to a structure residue.
    """
    def __init__(self,
                 column,
                 protein,
                 protein_res,
                 domain,
                 domain_res,
                 structure,
                 structure_res,
                 structure_dict):
        self.column = column
        self.protein = protein
        self.protein_res = protein_res
        self.domain = domain
        self.domain_res = domain_res
        self.structure = structure
        self.structure_res = structure_res

        # Make a copy. Don't keep references to the structure object.
        # It's potentially heavy weight.
        self.structure_dict = {}
        self.structure_dict.update(structure_dict)

    def index(self, index):
        """
        Visitor patter, the attribute will index itself.
        """
        index.column_protein[self.column].add((self.protein,self.protein_res))
        index.protein_domain[(self.protein.id,self.protein_res)] = (self.domain,self.domain_res)
        index.domain_structure[(self.domain.id,self.domain_res)].add((self.structure,self.structure_res))
        index.structure[(self.structure.index, self.structure_res)] = self

class IeaIndex(object):

    def __init__(self, family):
        self.mf = 0
        self.bp = 0
        self.cc = 0
        self.protein_term = defaultdict(lambda: set())
        self.term_protein = defaultdict(lambda: set())
        for protein in family.proteins:
            for iea in protein.iea:
                self.protein_term[protein].add(iea.term)
                self.term_protein[iea.term].add(protein)
                if iea.term_type=="molecular_function":
                    self.mf = self.mf+1
                if iea.term_type=="biological_process":
                    self.bp = self.bp+1
                if iea.term_type=="cellular_component":
                    self.cc = self.cc+1

    def terms(self):
        return self.term_protein.keys()

    def proteins(self, term):
        return self.term_protein[term]

class ColumnStructureIndex(object):

    def __init__(self):
        self.column_protein = defaultdict(lambda: set())
        self.protein_domain = {}
        self.domain_structure = defaultdict(lambda: set())
        self.structure = {}        
        self.dssp = defaultdict(lambda: None)
    
    def add(self, attribute):
        attribute.index(self)
    def columns(self):
        return sorted(set(self.column_protein.keys()))
    def proteins(self, column):
        return self.column_protein[column]
    def domain(self, protein, res):
        return self.protein_domain[(protein.id,res)]
    def structures(self, domain, res):
        return self.domain_structure[(domain.id,res)]
    def attributes(self, structure, res):
        return self.structure[(structure.index,res)]
    def getDssp(self,pdbid):
        return self.dssp[pdbid]
    def setDssp(self,pdbid,dssp):
        self.dssp[pdbid]=dssp

class FamilyFeatureBuilder(object):
    """
    Directs the provision of events for generating Family and structure feature documents.
    """

    def __init__(self,
                 handlerFactory,
                 structureProviderFactory,
                 columnProviderFactory,
                 ieaProviderFactory,
                 selectionProviderFactory):
        """
        @param handlerFactory: Factory method or classname for an xml.sax.ContentHandler.
        """
        self.handler = handlerFactory()
        self.structureProvider = structureProviderFactory(self.handler)
        self.columnProvider = columnProviderFactory(self.handler)
        self.ieaProvider = ieaProviderFactory(self.handler)
        self.selectionProvider = selectionProviderFactory(self.handler)

    def buildDocument(self, family):
        self.handler.startElement("family",Attributes(id=family.id))

        # Name, number of columns, etc
        self.handler.startElement("name")
        self.handler.characters(family.name)
        self.handler.endElement("name")

        self.handler.startElement("length")
        self.handler.characters(family.alignment.alignment.get_alignment_length())
        self.handler.endElement("length")

        self.handler.startElement("members")
        self.handler.characters(len(family.proteins))
        self.handler.endElement("members")

        # Proteins and family stuff
        self.structureProvider.buildProteins(family)
        self.selectionProvider.buildSelection(family)
        self.ieaProvider.buildIea(family)

        # Columns and Structure
        self.columnProvider.buildColumns(self.structureProvider.index)
        self.handler.endElement("family")


class SelectionFeatureProvider(object):

    def __init__(self, handler):
        self.handler = handler

    def buildSelection(self, family):
        self.handler.startElement("selection")
        for model in family.codeml:
            self.buildModel(model)
        self.handler.endElement("selection")

    def buildModel(self, model):
        self.handler.startElement("model", Attributes(id=model.model))
        if len(model.positive_selection)>0:
            self.handler.startElement("sites")
            for ps in model.positive_selection:
                self.buildSite(ps)
            self.handler.endElement("sites")
        self.handler.endElement("model")

    def buildSite(self, site):
        self.handler.startElement("site")

        self.handler.startElement("column")
        self.handler.characters(site.column)
        self.handler.endElement("column")

        self.handler.startElement("probability")
        self.handler.characters(site.probability)
        self.handler.endElement("probability")

        self.handler.startElement("post_mean")
        self.handler.characters(site.post_mean)
        self.handler.endElement("post_mean")

        self.handler.startElement("stderr")
        self.handler.characters(site.stderr)
        self.handler.endElement("stderr")

        self.handler.endElement("site")

class IeaFeatureProvider(object):

    def __init__(self, handler):
        self.handler = handler

    def buildIea(self, family):
        self.index = IeaIndex(family)
        self.handler.startElement("ontology")
        self.handler.startElement("mf")
        self.handler.characters(self.index.mf)
        self.handler.endElement("mf")
        self.handler.startElement("bp")
        self.handler.characters(self.index.bp)
        self.handler.endElement("bp")
        self.handler.startElement("cc")
        self.handler.characters(self.index.cc)
        self.handler.endElement("cc")
        self.buildTerms()
        self.handler.endElement("ontology")

    def buildTerms(self):
        self.handler.startElement("terms")
        for term in self.index.terms():
            if term != None:
                self.buildTerm(term)
        self.handler.endElement("terms")

    def buildTerm(self, term):
        self.handler.startElement("term")

        self.handler.startElement("acc")
        self.handler.characters(term.acc)
        self.handler.endElement("acc")
        self.handler.startElement("name")
        self.handler.characters(term.name)
        self.handler.endElement("name")
        self.handler.startElement("term_type")
        self.handler.characters(term.term_type)
        self.handler.endElement("term_type")

        self.handler.startElement("proteins")
        for protein in self.index.proteins(term):
            self.handler.startElement("protein", Attributes(id=protein.id,
                                                            ac=protein.ac.ac))
            self.handler.endElement("protein")    
        self.handler.endElement("proteins")    
        self.handler.endElement("term")


class StructureFeatureProvider(object):
    """
    Generates Protein, Domain, Structure events for the content handler.
    Alternatively, generates the DSSP features and indexes column-structure
    attributes.
    """
    
    def __init__(self, handler):
        self.handler = handler
        # Keep a mapping of columns to structure mappers
        self.index = ColumnStructureIndex()

    def buildProteins(self, family):
        self.alignment = family.alignment.alignment
        self.handler.startElement("proteins")
        for protein in family.proteins:
            self.buildProtein(protein)
        self.handler.endElement("proteins")

    def buildProtein(self, protein):
        self.protein = protein
        self.handler.startElement("protein",Attributes(id=protein.id,
                                                       ac=protein.ac.ac))
        
        self.buildDomains(protein)
        self.handler.endElement("protein")

    def buildDomains(self, protein):
        for domain_num,domain in enumerate(protein.domains):
            self.buildDomain(domain,domain_num)
    
    def buildDomain(self, domain, domain_num):
        self.domain = domain
        self.handler.startElement("domain",Attributes(id=domain.id,
                                              num=domain_num,
                                              start=domain.region.start,
                                              stop=domain.region.stop))
        self.handler.startElement("domain_type")
        self.handler.characters(domain.domain_type)
        self.handler.endElement("domain_type")
        self.buildStructures(domain)
        self.handler.endElement("domain")

    def buildStructures(self, domain):
        self.handler.startElement("structures")
        if domain.domain_type in ('psiblast','fold_recognition'):
            self.buildPDB(domain)
        
        # Index the structures, we will report the mcm scores grouped by structure_key
        structure_mcmdata = defaultdict(lambda: list())
        for mcm in domain.mcmdata:
            structure_mcmdata[mcm.structure].append(mcm)

        for structure, mcmdata in structure_mcmdata.items():
            self.buildDecoy(structure, mcmdata)
        
        self.handler.endElement("structures")

    def buildDecoy(self, structure, mcmdata):
        self.handler.startElement("structure",Attributes(type="decoy"))
        self.handler.startElement("id")
        self.handler.characters(str(structure.id))
        self.handler.endElement("id")
        
        self.handler.startElement("scop")
        for mcm in mcmdata:
            self.buildSCCS(mcm.sccs,mcm.probability)
        self.handler.endElement("scop")

        try:
            self.buildDecoyIndex(structure)
        except Exception, e:
            print e
        self.handler.endElement("structure")

    def buildPDB(self, domain):
        if not domain.sccs:
            return
        self.handler.startElement("structure",Attributes(type="pdb"))
        self.handler.startElement("id")
        self.handler.characters(str(domain.sccs.pdbid))
        self.handler.endElement("id")
        self.handler.startElement("chain")
        self.handler.characters(str(domain.sccs.chain))
        self.handler.endElement("chain")


        # This isn't guarenteed to have a superfamily prediction
        if domain.sccs and domain.sccs.sccs != None:
            self.handler.startElement("scop")
            for sccs in domain.sccs.sccs.split(","):
                self.buildSCCS(sccs,1.0)
            self.handler.endElement("scop")
        try:
            self.buildPDBIndex(domain)
        except Exception, e:
            print e
        self.handler.endElement("structure")
 
    def buildSCCS(self, sccs, probability):
        self.handler.startElement("sccs", Attributes(id=sccs))
        #self.handler.startElement("name")
        #self.handler.endElement("name")
        self.handler.startElement("probability")
        self.handler.characters(str(probability))
        self.handler.endElement("probability")
        self.handler.endElement("sccs")


    def buildPDBIndex(self, domain):
        # Now we have to build the mapper and index the columns it covers
        # This maps protein residues to the PDB structure record.
        from sqlalchemy.orm import object_session
        from hpf.hddb.db import PDBSeqRes
        parent_id = int(domain.parent_id[3:])
        session = object_session(domain)
        pdbseqres = session.query(PDBSeqRes).filter(PDBSeqRes.sequence_key==parent_id).first()

        protein_atom_mapper = ProteinToPDBAtom(self.protein, domain, self.alignment, pdbseqres=pdbseqres)
        # This maps protein residues to domain residues
        protein_domain_mapper = ProteinDomainMapper(self.domain)
        # This maps protein residues to alignment column numbers
        protein_alignment_mapper = ProteinAlignmentMapper(self.protein, self.alignment)

        # Simple structure object to use for indices
        structure = Structure(protein_atom_mapper._pdbid[0:4],chain=protein_atom_mapper._pdbid[4])
        print >>sys.stderr, "Protein",structure.id,structure.chain,protein_atom_mapper._pdbid
        
        # Get the PDB object downloaded and parsed for the mapping.
        chain = get_crystal(protein_atom_mapper._pdbid, structure=True)

        # Run DSSP on the set if not already run.
        dssp = self.index.getDssp(protein_atom_mapper._pdbid)
        if dssp!=None:
            print >>sys.stderr, "Using dssp",protein_atom_mapper._pdbid
        else:
            print >>sys.stderr, "Running dssp",protein_atom_mapper._pdbid
        features = DSSPFeatureSet(protein_atom_mapper,
                                  chain,
                                  dssp = dssp)
        self.index.setDssp(protein_atom_mapper._pdbid, features.dssp)
        
        # Iterate over the DSSP features and index them.
        for protein_res, xtra in features:
            column = protein_alignment_mapper[protein_res]
            domain_res = protein_domain_mapper[protein_res]
            structure_res = protein_atom_mapper[protein_res]
            attribute = ColumnStructureAttribute(column,
                                                 self.protein,
                                                 protein_res,
                                                 self.domain,
                                                 domain_res,
                                                 structure,
                                                 structure_res,
                                                 xtra)
            # This will index all of the attributes for later lookup and serialization.
            self.index.add(attribute)
            


    def buildDecoyIndex(self, structure):
        # A Structure object that is easier to index and safe of references.
        s = Structure(structure.id)

        # This maps protein residues to alignment column numbers
        protein_alignment_mapper = ProteinAlignmentMapper(self.protein, self.alignment)

        # Simple mapper from Protein to Domain.  The decoy has a direct mapping to the
        # decoy residues.
        protein_domain_mapper = ProteinDomainMapper(self.domain)
        # The first and only chain has no id, indexed as a space (' ').
        # This is the Rosetta prediction
        chain = list(structure.pdb[0])[0]

        # The DSSP structure Features
        features = DSSPFeatureSet(protein_domain_mapper,
                                  chain)
        # Iterate over the DSSP features and index them.
        for protein_res, xtra in features:
            column = protein_alignment_mapper[protein_res]
            domain_res = protein_domain_mapper[protein_res]
            structure_res = domain_res
            attribute = ColumnStructureAttribute(column,
                                                 self.protein,
                                                 protein_res,
                                                 self.domain,
                                                 domain_res,
                                                 s,
                                                 structure_res,
                                                 xtra)
            # This will index all of the attributes for later lookup and serialization.
            self.index.add(attribute)


class ColumnFeatureProvider(object):
    
    def __init__(self, handler):
        """
        """
        self.handler = handler
        
    def buildColumns(self,columnStructureIndex):
        self.index = columnStructureIndex
        self.handler.startElement("columns")
        for column in self.index.columns():
            self.buildColumn(column)
        self.handler.endElement("columns")

    def buildColumn(self, column):
        self.column = column
        self.handler.startElement("column", Attributes(id=column))
        self.buildProteins(column)
        self.handler.endElement("column")

    def buildProteins(self, column):
        self.handler.startElement("proteins",Attributes(id=column))
        for protein,protein_res in self.index.proteins(column):        
            self.buildProtein(protein,protein_res)
        self.handler.endElement("proteins")

    def buildProtein(self, protein, protein_res):
        self.handler.startElement("protein",Attributes(id=protein.id,
                                                       ac=protein.ac.ac,
                                                       residue=protein_res))
        domain, domain_res = self.index.domain(protein,protein_res)
        self.buildDomain(domain, domain_res)
        self.handler.endElement("protein")
        
    def buildDomain(self, domain, domain_res):
        self.handler.startElement("domain",Attributes(id=domain.id,
                                                      residue=domain_res))
        self.buildStructures(domain, domain_res)
        self.handler.endElement("domain")

    def buildStructures(self, domain, domain_res):
        self.handler.startElement("structures")
        for structure, structure_res in self.index.structures(domain,domain_res):
            self.buildStructure(structure, structure_res)
        self.handler.endElement("structures")

    def buildStructure(self, structure, structure_res):
        self.handler.startElement("structure",Attributes(residue=structure_res,
                                                         type=structure.type,
                                                         id=structure.index))
        self.handler.startElement("id")
        self.handler.characters(structure.id)
        self.handler.endElement("id")
        if structure.chain:
            self.handler.startElement("chain")
            self.handler.characters(structure.chain)
            self.handler.endElement("chain")

        att = self.index.attributes(structure, structure_res)
        for key, value in att.structure_dict.items():
            self.handler.startElement(key)
            self.handler.characters(value)
            self.handler.endElement(key)

        self.handler.endElement("structure")
