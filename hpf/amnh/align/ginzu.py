'''
Mappings for PSI-Blast and FFAS ginzu domains with AMNH families.
'''

from hpf.amnh.align import SeedAlignmentFactory, SeedAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align.Generic import Alignment

"""
Use a seed alignment to produce a domain to PDB mapping according to ginzu
results.
"""
class GinzuAlignmentFactory(SeedAlignmentFactory):
    
    def __init__(self):
        super(GinzuAlignmentFactory, self).__init__(_class = GinzuSeedAlignment)


    def _pdb(self,pdb_chain, domain_region):
        #TODO verify fencepost issues
        # The subsection of the PDB chain that ginzu reported an alignment to.
        pdb_start = domain_region.parent_start
        pdb_stop = domain_region.parent_stop

        # The pdb id looks like "P1;1m2vB/1-169"
        # We update the residue numbers to be the subsection we took.
        idparts = pdb_chain.id.split("/")
        pdb_id = "/".join([idparts[0],"%i-%i" % (pdb_start, pdb_stop)])
        pdb_section = SeqRecord(Seq(str(pdb_chain.seq)[pdb_start:pdb_stop],
                                    alphabet=pdb_chain.seq.alphabet), 
                                id=pdb_id, 
                                name=pdb_chain.name,
                                description=pdb_chain.description
                                )
        return pdb_section

    def _domain_alignment(self,alignment,domain_region, alignment_index):
        # Now we need to subselect the portion of the alignment 
        # that contains the domain.
        protein_record = alignment[alignment_index]
        protein_seq = str(protein_record.seq)
        # Figure out which columns encapsulate the domain.
        aa_count = 0
        column_start = None
        column_stop = None
        #print protein_seq
        for column,aa in enumerate(protein_seq):
            #print column,aa
            if aa!='-':
                aa_count=aa_count+1
            if aa_count==domain_region.start and column_start==None:
                column_start = column
            if aa_count==domain_region.stop and column_stop==None:
                column_stop = column
                break
        #print column_start,column_stop
        assert column_start != None, str(column_start)
        assert column_stop != None, str(column_stop)
        domain_alignment = Alignment(alphabet = alignment._alphabet)
        # Grab the portion of each sequence that correspond to columns
        # for the domain.
        for record in alignment:
            domain_alignment.add_sequence(record.id,
                                          str(record.seq)[column_start:column_stop])
        return (domain_alignment, column_start, column_stop)
        
    
    def create(self, alignment, domain_region, alignment_index, pdb_chain, **kwargs):
        """
        @param alignment: The seed alignment to use.
        @param domain_region: The ginzu domain information.
        @type domain_region: hpf.hddb.db.DomainRegion
        @param alignment_index: The sequence index in the alignment that corresponds to this target domain's entry.
        @param pdb_chain: The PDB Chain to align with.
        @type pdb_chain: Bio.SeqRecord.SeqRecord
        @return: SeedAlignment
        """
        pdb_section = self._pdb(pdb_chain,domain_region)
        domain_alignment,column_start, column_stop = self._domain_alignment(alignment,
                                                                            domain_region,
                                                                            alignment_index)
        ginzu_alignment = super(GinzuAlignmentFactory,self).create(domain_alignment,[pdb_section],**kwargs)
        ginzu_alignment.column_start = column_start
        ginzu_alignment.column_stop = column_stop
        ginzu_alignment.pdb_chain = pdb_chain
        ginzu_alignment.domain_region = domain_region
        ginzu_alignment.alignment_index = alignment_index
        return ginzu_alignment
        
                
class GinzuSeedAlignment(SeedAlignment):

    def __init__(self, *args, **kwargs):
        super(GinzuSeedAlignment,self).__init__(*args,**kwargs)
            
    def map(self,target_num=0):
        """
        
        """
        mapping = super(GinzuSeedAlignment,self).map(target_num)
        def column(m):
            return m+self.column_start if m!=None else None
        def nones(x):
            return [None for i in xrange(x)]
        return nones(self.domain_region.parent_start)+map(column,mapping)
        

if __name__=="__main__":
    from hpf.hddb.db import *
    session = Session()
    family = session.query(Family).get(307)
    al = family.alignments[0].alignment
    domains = list(session.query(Domain).filter(Domain.parent_sequence_key==304907).all())
    from hpf.pdb import *
    pdbid = '1yx6b'
    pdb = get_crystal(pdbid)
    pdb_seq = get_seq(pdb)
    from hpf.amnh.align import ginzu
    factory = ginzu.GinzuAlignmentFactory()
    alignment = factory.create(al,domains[0].region,0,SeqRecord(pdb_seq,pdbid))
    print str(alignment)
    print alignment.map()
    
