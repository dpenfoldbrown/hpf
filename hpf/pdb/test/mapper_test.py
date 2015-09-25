from hpf.hddb.db import *
from hpf.pdb.mapper import *
Session = Session()

family_id = 10187
domain_id = 753074

protein_id = 1060980
protein = Session.query(Protein).get(protein_id)

#kdrew: pulled this example from the full family alignment lib
def alignment_sample():

	family = Session.query(Family).get(family_id)
	al = family.alignment.alignment

	domain = Session.query(Domain).get(domain_id)



	# mapper = PDBDomainAlignmentMapper(protein,domain,al)
	# print "SEED TARGETS",mapper._seed._targets
	# list(mapper.mapping())
	# print "\n\n\n\n\n\n\n\n"
	# print "pdbmapper alignment"
	# for i in mapper.alignment():
	# 	print i.id
	# 	print str(i.seq)
	# print "\n\n\n\n\n\n\n\n"

	# mapper2 = PDBAtomMapper(protein,domain,mapper.alignment())
	# print "\n\n\n\n\n\n\n\n"
	# print "pdbatommapper alignment"
	# for i in mapper2.alignment():
	# 	print i.id
	# 	print str(i.seq)
	# print "\n\n\n\n\n\n\n\n"


	# print "Domain", mapper2._seed._records[mapper2._alignment_index].id
	# print str(mapper2._seed._records[mapper2._alignment_index].seq)
	# print "SequenceRecord", mapper2._seed._records[-2].id
	# print str(mapper2._seed._records[-2].seq)
	# print "AtomRecord", mapper2._seed._records[-1].id
	# print str(mapper2._seed._records[-1].seq)

	# for index_map in mapper2.mapping():
	# 	print index_map

        mapper = ProteinToPDBAtom(protein, domain, family.alignment.alignment)
	print 1,mapper._dict_()
	print 2,list(mapper.protein_to_alignment)[140:150]
	print 3,list(mapper.alignment_to_seqres)[0:10]
	print 4,list(mapper.seqres_to_atom)[0:10]

	print 176
	al = mapper.protein_to_alignment[176]
	print al
	s = mapper.alignment_to_seqres[al]
	print s
	a = mapper.seqres_to_atom[s]
	print a

if __name__=="__main__":
	alignment_sample()

