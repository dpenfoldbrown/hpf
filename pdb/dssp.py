from Bio.PDB.AbstractPropertyMap import AbstractResiduePropertyMap 
from Bio.PDB.DSSP import MAX_ACC
from Bio.PDB import to_one_letter_code, Structure, Model, Chain
import os
import re
_dssp_cys=re.compile('[a-z]')

SS = "SS_DSSP"
ASA = "EXP_DSSP_ASA"
RASA = "EXP_DSSP_RASA"

DESCRIPTION = {SS:"Secondary structure as calculated by DSSP",
               ASA:"Solvent accessibility calculated by DSSP",
               RASA:"Relative solvent accessibility calculated by DSSP"}

SHORT_NAME = {"rasa":RASA,
              "asa":ASA,
              "ss":SS}

def make_dssp_dict(input):
    """
    Return a DSSP dictionary that maps (chainid, resid) to
    aa, ss and accessibility, from DSSP output.

    @param input: the DSSP output
    @type input: string, open file handle
    """
    dssp={}
    if isinstance(input,basestring):
        fp= open(filename, "r")
    else:
        fp = input
    start=0
    keys=[]
    for l in fp.readlines():
        sl=l.split()
        if sl[1]=="RESIDUE":
            # start
            start=1
            continue
        if not start:
            continue
        if l[9]==" ":
            # skip -- missing residue
            continue
        resseq=int(l[5:10])
        icode=l[10]
        chainid=l[11]
        aa=l[13]
        ss=l[16]
        if ss==" ":
            ss="-"
        acc=int(l[34:38])
        res_id=(" ", resseq, icode)
        dssp[(chainid, res_id)]=(aa, ss, acc)
        keys.append((chainid, res_id))
    # close the file if we opened it
    if isinstance(input,basestring):
        fp.close()
    return dssp, keys


def dssp_dict_from_structure(structure, dssp="dsspcmbi", id=None):
    """
    In keeping with the spirit of functional programming...
    Creates a DSSP dictionary from more than just a filename.
    Pipes the structure into a dssp process if you send it one.
    @type structure: filename (string), Bio.PDB.Entity, hpf.hddb.db.Structure
    @return (Bio.PDB.Entity, dssp, keys)
    """
    from subprocess import Popen,PIPE
    from cStringIO import StringIO

    dssp = Popen("which %s"%dssp, shell=True, stdout=PIPE).communicate()[0].strip()
    assert os.path.exists(dssp)
    # If we have a structure object we will pipe it into DSSP
    # otherwise it will be loaded and parsed.
    if isinstance(structure,basestring):
        pipe=None
        dssp = "%s %s" % (dssp, structure)
        from Bio.PDB import PDBParser
        structure = PDBParser().get_structure(id if id else structure,structure)       
    else:
        # We write to a buffer which will be piped into the DSSP stdin
        # using pipe.getvalue(), the dssp stdin option must be added
        dssp = "%s --" % dssp
        import tempfile
        pipe = tempfile.NamedTemporaryFile('w')
        from Bio.PDB.Structure import Entity, Structure as BioStructure
        if isinstance(structure,Entity):
            s = structure
            while(s.get_parent()!=None):
                s = s.get_parent()
            from Bio.PDB.PDBIO import PDBIO
            io = PDBIO()
            io.set_structure(s)
            io.save(pipe.name)
        else:
            # Otherwise this is of type hpf.hddb.db.Structure
            # avoid importing the class and mapped meta table stuff
            with pipe as handle:
                pipe.write(structure.text)
            structure = structure.structure

    # execute the process, piping data if necessary
    #print pipe.getvalue()
    out,err= Popen(dssp,shell=True,stdin=pipe if pipe else None,stdout=PIPE,stderr=PIPE).communicate()
    io = StringIO()
    io.write(out)
    io.seek(0)
    return (structure,)+make_dssp_dict(io)

class DSSP(AbstractResiduePropertyMap):
    """
    Modifies the DSSP annotator class to catch exceptions and skip residues with problems.
    Fails silently on residues.  Accepts hpf.hddb.db.Structure objects that can be piped into DSSP.
    @see Bio.PDB.DSSP.DSSP
    """
    def __init__(self, structure, dssp="dsspcmbi", id=None):
        """
        @param model: the first model of the structure
        @type model: L{Model}

        @param pdb_file: a PDB file
        @type pdb_file: string

        @param dssp: the dssp executable (ie. the argument to os.system)
        @type dssp: string

        @param id: the structure id to use if parsed
        @type id: string
        """
        # create DSSP dictionary
        structure, dssp_dict, dssp_keys = dssp_dict_from_structure(structure, dssp, id=id)
        self.structure = structure

        self.model = structure
        while self.model.get_parent() != None:
            self.model = self.model.get_parent()
        self.model = self.model[0]
        #assert isinstance(self.model, Model)

        dssp_map={}
        dssp_list=[]
        # Now create a dictionary that maps Residue objects to 
        # secondary structure and accessibility, and a list of 
        # (residue, (secondary structure, accessibility)) tuples
        for key in dssp_keys:
            try:
                chain_id, res_id=key
                
                chain=self.model[chain_id]
                res=chain[res_id]
            except:
                continue
            aa, ss, acc=dssp_dict[key]
            res.xtra[SS]=ss
            res.xtra[ASA]=acc
            # relative accessibility
            resname=res.get_resname()
            rel_acc=acc/MAX_ACC[resname]
            if rel_acc>1.0:
                rel_acc=1.0
            res.xtra[RASA]=rel_acc
            # Verify if AA in DSSP == AA in Structure
            # Something went wrong if this is not true!
            resname=to_one_letter_code[resname]
            if resname=="C":
                # DSSP renames C in C-bridges to a,b,c,d,...
                # - we rename it back to 'C'
                if _dssp_cys.match(aa):
                    aa='C'
            if not (resname==aa):
                raise PDBException("Structure/DSSP mismatch at "+str(res))
            dssp_map[key]=((res, ss, acc, rel_acc))
            dssp_list.append((res, ss, acc, rel_acc))
        AbstractResiduePropertyMap.__init__(self, dssp_map, dssp_keys, dssp_list)

    def get_structure(self):
        """
        @return The Bio.PDB.Entity passed as the init structure, or the 
        Entity loaded from the init structure object.
        """
        return self.structure

    def get_model(self):
        """
        @return The first model from the structure supplied.
        """
        return self.model

class DSSPFeatureSet(object):
    """
    Calculate DSSP features based on some domain structure and map them to Protein indices.
    """
    def __init__(self,
                 proteinStructureMapper,
                 chain,
                 dssp=None):
        """
        @param proteinStructureMapper: A mapping from protein indices to structure indices.
        @type proteinStructureMapper: hpf.mapper.MapperInterface
        @param chain: A Biopython Chain object.  Does not refer to a PDB chain, per se, but is
        just part of the biopython structure object heirarchy.  Even Rosetta models, loaded in biopython
        will have one chain.
        @type chain: Bio.PDB.Chain
        @param key: The type of DSSP info to return.
        @param dssp: A DSSP object, if provided DSSP will not be run on the structure.
        """
        self.mapper = proteinStructureMapper
        self.chain = chain
        # Annotate the structure with DSSP features
        self.dssp = DSSP(chain) if dssp==None else dssp

    def _key(self,key):
        if key not in (SHORT_NAME.values()):
            key = SHORT_NAME[key]
        return key

    def description(self, key):
        return DESCRIPTION[self._key(key)]

    def scores(self, key):
        for protein,xtra in self:
            if xtra.has_key(self._key(key)):
                yield (protein,xtra[self._key(key)])

    def __iter__(self):
        """
        Yields protein indices paired with DSSP scores.
        @return generator of (protein,chain[protein].xtra[key]) tuples.
        """
        atom_res = list(self.chain)
        for protein,structure in self.mapper:
            if protein==None or structure==None:
                continue
            xtra = atom_res[structure].xtra
            yield (protein, xtra)

if __name__=="__main__":
    pdbid="1n6j"
    from hpf.pdb import get_crystal
    pdb = get_crystal("1yelA", structure=True, pdir="/tmp")
    dssp = DSSP(pdb)
    print "done"
    print pdb
    for res in pdb:
        rasa = res.xtra[RASA] if res.xtra.has_key(RASA) else None
        print res.get_resname()[1], rasa

    from hpf.hddb.db import Session, Structure
    session = Session()
    pdb = list(session.query(Structure).get(1054611).pdb[0])[0]
    dssp = DSSP(pdb)
    for res in pdb:
        rasa = res.xtra[RASA] if res.xtra.has_key(RASA) else None
        print res.get_resname()[1], rasa
