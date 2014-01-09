import os
import sys
import subprocess
from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO

class AstralList(object):
     
    def __init__(self,
                 pdir='/scratch/kdrew/data/astral/models', 
                 sdir='/scratch/kdrew/data/astral/fasta', 
                 obsolete_pdb=None):
        self.pdir = pdir
        self.sdir = sdir
      
    def retrieve_astral_file(self, astral_id, obsolete=0, pdir=None):
        """Retrieves a PDB structure file from the PDB server and stores it in a local file tree"""
        if not pdir:
            pdir = self.pdir
        file = os.path.join(pdir,"%s.ent" % astral_id)

        if os.path.exists(file) and os.path.getsize(file) > 0:
            pass
        else:
            url = "http://astral.berkeley.edu/pdbstyle.cgi?id=%s&output=text" % astral_id
            cmd = "wget -O %s \"%s\"" % (file,url)
            print >>sys.stderr, "Astral: d/l %s" % astral_id
            with open(os.devnull,"w") as null:
                subprocess.check_call(cmd,shell=True,stdout=null,stderr=null)
        return file

    def retrieve_astral(self, astral_id, pdir=None):
        file = self.retrieve_astral_file(astral_id,pdir=pdir)
        return PDBParser().get_structure(astral_id,file)
    
    def retrieve_astral_sequence_fasta(self,astral_id,sdir=None,atom=False):
        """Retrieves an astral fasta file from the server and stores it in a local file tree"""
        if not sdir:
            sdir = self.sdir
        file = os.path.join(sdir,"%s_%s.fasta" % (astral_id, "atom" if atom else "seqres"))

        if os.path.exists(file) and os.path.getsize(file) > 0:
            pass
        else:
            url = "http://astral.berkeley.edu/getseqs/getseqs.cgi?id=%s&kind=GD&seqOption=%i&output=text" % (astral_id, 1 if atom else 0)
            cmd = "wget -O %s \"%s\"" % (file,url)
            print >>sys.stderr, "Astral Seq: d/l %s" % astral_id
            print cmd
            with open(os.devnull,"w") as null:
                subprocess.check_call(cmd,shell=True,stdout=null,stderr=null)
        return file
    
    def retrieve_astral_sequence_record(self, astral_id, sdir=None, atom=0):
        file = self.retrieve_astral_sequence_fasta(astral_id,sdir=sdir,atom=atom)
        with open(file) as handle:
            return list(SeqIO.parse(handle, "fasta"))[0]

#download_entire_pdb(self, listfile=None)
#Retrieves all PDB entries not present in the local PDB copy.     source code
#      
#download_obsolete_entries(self, listfile=None)
#Retrieves all obsolete PDB entries not present in the local obsolete PDB copy.     source code
#      
#get_all_entries(self)
#Retrieves a big file containing all the PDB entries and some annotation to them.     source code
#      
#get_all_obsolete(self)
#Returns a list of all obsolete entries ever in the PDB.     source code
#      
#get_recent_changes(self)
#Returns three lists of the newest weekly files (added,mod,obsolete).     source code
#      
#get_seqres_file(self, savefile='pdb_seqres.txt')
#Retrieves a (big) file containing all the sequences of PDB entries and writes it to a file.     source code
#      
#get_status_list(self, url)
#Retrieves a list of pdb codes in the weekly pdb status file from the given URL.     source code
#      
#      
#update_pdb(self)
#I guess this is the 'most wanted' function from this module.
