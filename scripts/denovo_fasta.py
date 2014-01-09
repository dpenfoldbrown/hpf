
from hpf_report.model.meta import Session
from hpf.hddb.db import Experiment, Domain, Protein, Sequence, DomainSCCS, YRC
from itertools import imap
from sqlalchemy import not_
from sqlalchemy.orm.exc import NoResultFound


from Bio import SeqIO

class YRCRecordFactory(object):
	"""Create a Bio.SeqRecord from an hddb Sequence object"""
	def create(self, *sequences, **kwargs):
		from Bio.SeqRecord import SeqRecord
		from Bio.Seq import Seq
		for sequence in sequences:
			yrc = Session().query(YRC).filter(YRC.sequence_key == sequence.id).first()
			if yrc:
				description="yrc_sequence_key:"+str(yrc.yrc_sequence_key)+"|yrc_protein_key:"+str(yrc.yrc_protein_key)
				#description+=",".join([y.yrc_protein_key for y in yrc])
				print description
			else:
				print "no yrc info: ", sequence.id
				description = ""
			yield SeqRecord(Seq(sequence.sequence),str(sequence.id), description=description,**kwargs)


def domain_fasta(id=None):

	if id: 
		outfile = "/Users/kdrew/tmp/exp_%d.fasta" % (id,)
		experiment = Session.query(Experiment).filter(Experiment.id==id).one()
		species = experiment.species()
		print species
		sequences = Session().query(Sequence).join(DomainSCCS.domain).join(Domain.sequence).join(Domain.proteins).join(Protein.experiment).filter(not_(DomainSCCS.domain_type.in_( ('psiblast','fold_recognition')))).filter(Experiment.id==id).distinct().all()
	else:
		print "all denovo"
		outfile = "/Users/kdrew/tmp/hpf_denovo.fasta"
		sequences = Session().query(Sequence).join(DomainSCCS.domain).join(Domain.sequence).join(Domain.proteins).join(Protein.experiment).filter(not_(DomainSCCS.domain_type.in_( ('psiblast','fold_recognition')))).distinct().all()


	print len(sequences)

	#sequences = Session().query(Sequence).join(Domain.sequence).join(Domain.proteins).join(Protein.experiment).filter(Experiment.id==id).distinct().all()
	#records = imap(lambda x: x.biopython(description=species).format("fasta"), sequences)


	bio_records = YRCRecordFactory().create(*sequences)
	with open(outfile, "w") as handle1:
		SeqIO.write(bio_records, handle1,"fasta")

	handle1.close()
        
#domain_fasta(id=  1163)
#domain_fasta(id= 805)
domain_fasta()

