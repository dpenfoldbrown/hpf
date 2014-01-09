"""
@author kdrew
Mapper module for mapping alignments to pdb sequences including atom records
"""
from hpf.mapper import DomainAlignmentMapper, MapperInterface, MapperChain, ProteinAlignmentMapper
from hpf.hddb.db import Session, Family, Protein, Domain, PDBSeqRes

class AlignmentToPDBMapper():
# Creates a mapper between a family alignment and a pdb structure
# NOTE: This is an abstraction of all of the confused and undocumented mappers below.
# Use this, because it is clear and because it works. KEEP IT SIMPLE.
# UNIT TESTS in pdb/tests/alignment_pdbmapper_test.py
# Access three public maps:
#   alignment_pdbseq_map  - alignment col => pdb seqres
#   pdbseq_pdbatom_map    - pdb seqres    => pdb atom
#   alignment_pdbatom_map - alignment col => pdb atom

    def __init__(self, family, protein, domain, debug=False):
        self.DEBUG = debug
        self.session = Session()
        self.family = family
        self.protein = protein
        self.domain  = domain
        self.alignment = self.family.alignments[0]
        if self.domain.sccs:
            self.pdbseqres = self.session.query(PDBSeqRes).filter_by(sequence_key=self.domain.parent_id[3:], chain=self.domain.sccs.chain).first()
        else:
            self.pdbseqres = self.session.query(PDBSeqRes).filter_by(sequence_key=self.domain.parent_id[3:]).first()
        self.pdbid = self.pdbseqres.pdb.pdbId+self.pdbseqres.chain

        # Create map: alignment column -> pdb seq res
        self.alignment_pdbseq_map = PDBDomainAlignmentMapper(self.protein, self.domain, self.alignment.alignment, self.pdbseqres, inverse=True)

        # Create map: pdb seq res -> pdb atom res
        self.pdbseq_pdbatom_map = PDBAtomSeqResMapper(self.pdbseqres, inverse=False)

        # Create map: alignment column -> pdb atom res
        self.alignment_pdbatom_map = AlignmentToPDBAtomMapper(self.alignment_pdbseq_map, self.pdbseq_pdbatom_map)


class PDBDomainAlignmentMapper(DomainAlignmentMapper):
# MAPPING: PDB Seqres -> Family alignment column (inverse gives alignment column -> pdb seqres)
    """
    Map a domain's PDB GINZU match using it's database Sequence Record to 
    a family constraining the alignment to land on the target domain.
    This is higher resolution than mapping only to domain residues, since the
    PDB is actually aligned to the family's alignment.  The alignment may be
    gapped in the domain's position, but the family could still have selection
    analysis used to paint the structure.
    """
    
    def __init__(self, 
             protein, 
             domain, 
             alignment, 
             pdbseqres=None, 
             **kwargs):
        """
        @type protein: hpf.hddb.db.Protein
        @type domain: hpf.hddb.db.Domain
        @type alignment: Bio.Align.Generic
        @param pdbseqres: (Optional) Specify a PDB chain's sequence.
        """
        super(PDBDomainAlignmentMapper,self).__init__( protein, domain, alignment, **kwargs)
        self._chain = pdbseqres
        
        # If no pdbseqres given, fetch PDBSeqRes ORM object from db by sequence_key (a parsed domain.parent_id). 
        if self._chain==None:
            from hpf.hddb.db import Session,PDBSeqRes
            Session = Session()
            parent_id = int(self._domain.parent_id[3:])
            if domain.sccs:
                self._chain = Session.query(PDBSeqRes).filter_by(sequence_key=parent_id, chain=domain.sccs.chain).first()
            else:
                self._chain = Session.query(PDBSeqRes).filter(PDBSeqRes.sequence_key==parent_id).first()
        pdbid = self._chain.pdb.pdbId+self._chain.chain
        self._pdbid = pdbid;
        self._seed = self._seed_alignment()

    def get_structure(self):
        """
        @return Bio.PDB.Structure from the PDB's chain
        """
        return self._chain.structure.structure

    def get_seqres(self):
        return self._chain

    def _seed_alignment(self):
        """
        Generate a seed alignment that includes the 
        PDB's sequence record.
        """
        # Build a Seed Alignment
        from hpf.amnh.align import ginzu
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC

        from Bio.Seq import Seq

        #parent_id = int(self._domain.parent_id[3:])
        # use the sequence given in the chain object
        pdb_sequence = self._chain.sequence
        #print "pdb_sequence:", pdb_sequence.sequence

        factory = ginzu.GinzuAlignmentFactory()
        # The factory creates a subalignment that only
        # includes columns within the domain region
        alignment = factory.create(self._alignment,
                       self._domain.region,
                       self._alignment_index,
                       pdb_sequence.record)
        #SeqRecord(Seq(pdb_sequence.sequence,IUPAC.protein),self._pdbid))
        return alignment

    def alignment(self):
        """
        Returns an alignment showing the domain aligned with the PDB.
        """
        pdb = self._seed._records[-1]
        #print pdb
        domain = self._seed._records[self._alignment_index]
        # Remove columns gapped in both.  These are gapped because
        # the seed alignment includes the entire family.
        pdb_seq = []
        domain_seq = []
        for p,d in zip(str(pdb.seq),str(domain.seq)):
            if(p=="-" and d=="-"):
                continue
            pdb_seq.append(p)
            domain_seq.append(d)
        from Bio.Align.Generic import Alignment
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        start = str(self._domain.region.start)
        stop = str(self._domain.region.stop)
        align = Alignment(self._seed._alphabet)
        #name = self._protein.ac.ac+"_"+start+":"+stop
        name = str(self._protein.id)
        description = "Domain identified between protein residues "+start+" and "+stop
        align._records.append(SeqRecord(Seq("".join(domain_seq)),id=name,description=description))
        align._records.append(SeqRecord(Seq("".join(pdb_seq)),id=pdb.id))
        return align

    def mapping(self):
        """
        Yields a (residue,column) tuple generator that maps PDB
        residues to columns in the alignment.
        """
        for residue,mapping in enumerate(self._seed.map()):
            if mapping!=None:
            # NOTE: ADDING THE +1 to enumerated residue number because PDB Seqres records start from 1, not 0
            # This yields a correct mapping (I think. Check with own results, but for +Sel this fits)
                yield (residue+1, mapping)

class PDBAtomMapper(PDBDomainAlignmentMapper):
    """
    @deprecated 
    dpb: Unsure of what mapping this returns. Is [something] -> alignment column.
    """
    def _seed_alignment(self):
        # Build a Seed Alignment
        from hpf.pdb import get_seq, get_crystal
        #kdrew: this gets the atom record sequence and not the full sequence
        pdb_seq = get_seq(get_crystal(self._pdbid))
        #print "_seed_alignment pdb sequence: ", pdb_seq
        from hpf.amnh.align import SeedAlignmentFactory
        from Bio.SeqRecord import SeqRecord
        #kdrew: might change this to be pairwise2 alignment
        factory = SeedAlignmentFactory()
        alignment = factory.create(self._alignment,
                       [SeqRecord(pdb_seq,self._pdbid+"_atom_record")])
        return alignment


class AlignmentToPDBAtomMapper(MapperInterface):
# A family alignment to pdb atom record mapper.
# Takes:
#   PDBDomainAlignmentMapper (inverse=True) al_to_pdbseq    - an alignment column -> pdb seqres map
#   PDBAtomSeqResMapper (inverse=False) pdbseq_to-pdbatom   - a pdb seqres -> pdb atom map
# MAPPING: family alignment column -> pdb atom record

    def __init__(self, al_to_pdbseq, pdbseq_to_pdbatom, **kwargs):
        super(AlignmentToPDBAtomMapper,self).__init__(**kwargs)
        
        self.al_to_pdbseq = al_to_pdbseq
        self.pdbseq_to_pdbatom = pdbseq_to_pdbatom

    def mapping(self, ):
        for col, seqres in self.al_to_pdbseq:
            try:
                atom = self.pdbseq_to_pdbatom[seqres]
            except KeyError, e:
                #DEBUG print "PDB residue {0} does not exist in pdbseq to atom map".format(seqres)
                #DEBUG print e
                continue
            yield (col, atom)

    def get_structure(self, ):
        return self.pdbseq_to_pdbatom.get_structure()


class PDBAtomSeqResMapper(MapperInterface):
# MAPPING: PDB atom record -> PDB seqres (can invert)

    def __init__(self, pdbseqres, **kwargs):
        """
        @type chain: Bio.PDB.Chain
        """
        super(PDBAtomSeqResMapper,self).__init__(**kwargs)
        self.pdbseqres = pdbseqres

        # Parse PDBSeqRes resmap field (contains pdb seq -> pdb atom mapping). 
        self.pdbseq_to_pdbatom = self._parse_resmap(self.pdbseqres.resmap)


    def mapping(self):
        # Query to DB to get pdbseq->pdbatom map; parse map into list of tuples [(seq1, atom1), (seq2,atom2),...];
        # yield tuple vals (for seq, atom in list)
        # Takes a PDBSeqRes object
        # NOTE: This mapping does NOT contain pdb seqrecs that do not map to pdb atom recs. Will throw a key error 
        # if you attempt to access a seq->atom mapping that is not present.

        for pdbseq, pdbatom in self.pdbseq_to_pdbatom:
            if pdbatom == -9999:
                continue
            yield (pdbseq, pdbatom)

    def get_structure(self):
        return self.pdbseqres.structure

    def _parse_resmap(self, resmap):
        # Parses a resmap string of the form found in the hpf database, pdbSeqRes table, resmap column.
        # (According to regexp below)
        import re
        maplist = resmap.splitlines()
        map_pattern = r"(?P<res>[a-zA-Z]{1})\s+(?P<pdbseq>[0-9]+)\s+(?P<pdbatom>-?[0-9]+)"
        map_tuples = []
        for map in maplist:
            map_found = re.match(map_pattern, map)
            if not map_found:
                raise "Pattern '{0}' does not match row '{1}' of PDBSeqRes object {2}".format(map_pattern, map, self.pdbseqres.id)
            map_tuples.append((int(map_found.group('pdbseq')), int(map_found.group('pdbatom'))))
        return map_tuples


class ProteinToPDBAtom(MapperChain):
    """
    Creates a Chained Mapping from a protein to the PDB atom record
    for use in accessing structure objects.
    MAPPING: protein sequence -> PDB atom record
    """
    def __init__(self, protein, domain, alignment, pdbseqres=None, **kwargs):
        # Get pdbseqres if nothing given (needed for PDBAtomSeqResMapper)
        if pdbseqres == None:
            from hpf.hddb.db import Session, PDBSeqRes
            if domain.sccs:
                pdbseqres = Session().query(PDBSeqRes).filter_by(sequence_key=domain.parent_id[3:], chain=domain.sccs.chain).first()
            else:
                pdbseqres = Session().query(PDBSeqRes).filter_by(sequence_key=domain.parent_id[3:]).first()
        
        # Maps PDB Sequence Record -> Alignment
        self.alignment_to_seqres = PDBDomainAlignmentMapper(protein,domain,alignment,pdbseqres=pdbseqres, inverse=True)
        self._pdbid = self.alignment_to_seqres._pdbid

        # Maps Seqres to actual Atom record indices.
        self.seqres_to_atom = PDBAtomSeqResMapper(pdbseqres, inverse=True)

        # Maps Alignment -> Protein
        self.protein_to_alignment = ProteinAlignmentMapper(protein, alignment)
        self.alignment_to_protein = ProteinAlignmentMapper(protein, alignment, inverse=True)
        
        # Create the chain
        super(ProteinToPDBAtom,self).__init__(self.protein_to_alignment,
                              self.alignment_to_seqres,
                              self.seqres_to_atom)

    def get_structure(self):
        return self.seqres_to_atom.get_structure()

    def get_chain(self):
        return self.seqres_to_atom.get_chain()
