import Bio.PDB
from Bio.PDB.PDBIO import Select
from Bio.PDB.DSSP import make_dssp_dict,ss_to_index,MAX_ACC
import Bio.PDB.DSSP
import Bio.AlignIO.ClustalIO

import subprocess
import os
import tempfile

DSSP_DIR="/Users/patrick/experiment/data/dssp"

class RMSDSuperimposer(Bio.PDB.Superimposer):
    """Performs RMS analysis between two structures based on an structure alignment"""
    
    def __init__(self, alignment=None, structures=None):
        Bio.PDB.Superimposer.__init__(self)
        fixed, moving = ([],[])
        if alignment:
            for r1,r2 in alignment.get_iterator():
                print r1,r2
                fixed += list(r1)
                moving += list(r2)
        elif structures:
            s1,s2 = structures
            fixed += list(s1.get_atoms())
            moving += list(s2.get_atoms())
        self.set_atoms(fixed, moving)
        

class PDBAlignment(Bio.PDB.StructureAlignment):
    from Bio import SeqRecord
    """Performs a structure alignment based on the polypeptide sequences of two structure"""
    
    def __init__(self, structure1, structure2, ppb=Bio.PDB.PPBuilder() ):
        self.structure1 = structure1
        self.structure2 = structure2
        self.sequence1 = SeqRecord(get_seq(structure1, ppb=ppb),0)
        self.sequence2 = SeqRecord(get_seq(structure2, ppb=ppb),1)
        print self.sequence1, self.sequence2
        self.alignment = align_seq(self.sequence1, self.sequence2)
        Bio.PDB.StructureAlignment.__init__(self,self.alignment,self.structure1, self.structure2)

class PDBChainSelector(Select):
    
    def __init__(self, chain):
        self.chain = chain
        self.found = False
    
    def accept_chain(self, chain):
        r = chain.id == self.chain
        if r:
            self.found = True
        return r

class DSSP(dict):

    SS_PAIRS = ['HH','HB','HE','HG','HI','HT','HS','H-','BB','BE','BG','BI','BT','BS','B-','EE','EG','EI','ET','ES','E-','GG','GI','GT','GS','G-','II','IT','IS','I-','TT','TS','T-','SS','S-','--']
    #H        Alpha helix (4-12) 
    #B        Isolated beta-bridge residue 
    #E        Strand 
    #G        3-10 helix 
    #I        pi helix 
    #T        Turn 
    #S        Bend 
    #-        None 

    @staticmethod
    def index(s1,s2):
        """Return the secondary structure index for residues r1,r2"""
        if s1==None:
            s1='-'
        if s2==None:
            s2='-'
        ss = str(s1)+str(s2)
        try:
            return DSSP.SS_PAIRS.index(ss)
        except:
            try:
                return DSSP.SS_PAIRS.index(reverse(ss))
            except:
                raise Exception("NO SSPAIR",ss)

    def ss(self, residue):
        info = self.info(residue)
        if info:
            #print info
            aa,ss,acc = info
            return ss
        return None

    def rel_acc(self,residue):
        """Return the relative accessibility of the residue"""
        info = self.info(residue)
        if info:
            aa,ss,acc = info
            if acc==None:
                acc=0.0
            resname=residue.get_resname()
            rel_acc=acc/MAX_ACC[resname] 
            if rel_acc>1.0: 
                rel_acc=1.0 
            if rel_acc==None:
                rel_acc = 0.0
            return rel_acc
        return None

    def info(self, residue):
        key = residue.get_full_id()[-2:]
        if self.has_key(key):
            return self[key]
        else:
            return None

def get_dssp(pdbid):
    f = os.path.join(DSSP_DIR,"%s.dssp"%pdbid)
    assert(os.path.exists(f))
    dssp,keys = make_dssp_dict(f)
    return DSSP(dssp)

def get_crystal(id,pdir="/tmp/structures",structure=True):
    """
    Get a crystal structure, using the id to auto-determine if this is a pdbid,
        pdbid and chain, or astral id
    """
    if len(id)==4:
        return get_pdb(id, structure=structure, pdir=pdir)
    elif len(id)==5:
        return get_pdb(id[0:4], "A" if id[4]=="_" else id[4].upper(), structure=structure, pdir=pdir)
    else:
        return get_astral(id,structure=structure,pdir=pdir)

def get_pdb(pdbid, chain=None, structure=True, pdir="/scratch.tmp/pcw216/shared/pdblist"):
    """
    Returns a PDB structure downloaded from the rcsb.
    Uses the HTTP alternative since the FTP server denies too many connections
    from a single IP address, affecting parallel runs on clusters."""
    chain = chain.upper() if chain else chain
    chain = "A" if chain == "_" else chain
    pdbl = Bio.PDB.PDBList(server=Bio.PDB.PDBList.alternative_download_url,pdb=pdir)

    # Stupid thing prints to stdout, redirect temporarily
    import sys
    out = sys.stdout
    sys.stdout = sys.stderr
    f = pdbl.retrieve_pdb_file(pdbid,pdir=pdir)
    sys.stdout = out
    if structure:
        s = Bio.PDB.PDBParser().get_structure(pdbid,f)
        for c in s.get_chains():
            if chain != None and c.id==chain:
                return c
        if chain != None:
            raise Exception("Chain not found",chain)
        return s
    else:
        return f
    
def get_astral(astral_id, structure=True, pdir="/tmp/structures"):
    """Returns a PDB structure downloaded from the astral ftp server"""
    from hpf.pdb.astral import AstralList
    al =  AstralList(pdir=pdir, sdir=pdir, obsolete_pdb=True)
    if structure:
        return al.retrieve_astral(astral_id)
    else:
        return al.retrieve_astral_file(astral_id)
        
def align_seq(sequences, blast=False,mat='BLOSUM'):
    if blast:
        return align_seq_blast(sequences)
    else:
        return align_seq_clustalw(sequences)

def align_seq_blast(sequences):
    pass

def align_profile_seq(profile, fasta,mat='BLOSUM',output_order='INPUT',outfile=None):
    """Align sequences to a given profile."""
    if outfile:
        outfile = open(outfile,"w")
    else:
        outfile = tempfile.NamedTemporaryFile("w")
    try:
        cmd =  "clustalw -PROFILE1='%s' -PROFILE2='%s' -OUTFILE='%s' -OUTORDER='%s' -SEQUENCES" % (profile, fasta, outfile.name,output_order)
        subprocess.check_call(cmd,shell=True,stdout=open(os.devnull,"w"))
        f = open(outfile.name)
        try:
            alignments = list(Bio.AlignIO.ClustalIO.ClustalIterator(f))
        finally:
            f.close()
    finally:
        outfile.close()                                             
    return alignments

def align_seq_clustalw(sequences,mat='BLOSUM',output_order='INPUT'):
    """Return a clustal alignment of an arbitrary number of sequences.  Requires clustalw on path.
    Sequences are of type SeqRecord"""         
    assert(len(sequences) > 1)
    from Bio.Clustalw import MultipleAlignCL
    from Bio import Clustalw
    from Bio.SeqIO.FastaIO import FastaWriter

    import tempfile

    def write(f):
        i = 0
        for seqrecord in sequences:
            if seqrecord.id == "<unknown id>":
                seqrecord.id = str(i)
            i+=1
        writer = FastaWriter(f)
        writer.write_file(sequences)
        f.flush() #IMPORTANT
    
    input = tempfile.NamedTemporaryFile("w")
    #input = open("input","w")
    try:
        write(input)
        cline = MultipleAlignCL(input.name)
        cline.set_protein_matrix(MultipleAlignCL.PROTEIN_MATRIX[MultipleAlignCL.PROTEIN_MATRIX.index(mat)])
        cline.is_quick = False
        # Set to PAM matrix.
        output = tempfile.NamedTemporaryFile("w")
        #output = open("output","w")
        try:
            cline.set_output(output.name, output_order=MultipleAlignCL.OUTPUT_ORDER[MultipleAlignCL.OUTPUT_ORDER.index(output_order)])
            #print "input",open(input.name).read()
            #print "cline",cline
            subprocess.check_call(str(cline),shell=True,stdout=open(os.devnull,"w"))
            f = open(output.name)
            try:
                iter = Bio.AlignIO.ClustalIO.ClustalIterator(f)
                # There should only be one alignment in the file
                alignments = list(iter)
                assert(len(alignments)==1)
                alignment = alignments[0]
            finally:
                f.close()
            #print "output",output.read()
            return alignment
        finally:
            output.close()
    finally:
        input.close()


def get_pp(structure, ppb=None):
    """Extract the polypeptide from a structure, concatenating all chains"""
    if not ppb:
        ppb = Bio.PDB.CaPPBuilder()
    all = Bio.PDB.Polypeptide.Polypeptide()
    for pp in ppb.build_peptides(structure):
        all = all+pp
    return Bio.PDB.Polypeptide.Polypeptide(all)
    

def get_seq(structure, ppb=None):
    """Extract the sequence from a pdb structure.  It's highly suggested you only pass in single chains or you'll receive everything!
    @return: Bio.Seq.Seq
    """
    pp = get_pp(structure,ppb=ppb)
    return pp.get_sequence()

def res_dist(n, p):
    """Return the distance between the c-alpha's of two residues or c-beta"""
    return n['CA']-p['CA']

def reverse( thing ):
    """Reverse a list or tuple"""
    def revi( a, b ):
        if not a : 
            return b
        return revi( a[1:], a[:1]+b )
    return revi( thing, thing[:0] )
