import sys, os
from hpf.hddb.db import Experiment, Protein, Sequence, Session

session = Session()

def main():
	sequences = session.query(Sequence.id,Sequence.sequence).join(Protein.sequence, Protein.experiment).distinct()
	#sequences = session.query(Sequence.id,Sequence.sequence).join(Protein.sequence, Protein.experiment).filter(Protein.experiment_key == 901).distinct()
	for sequence in sequences:
		print ">bddb|%s| [BDDB-SEQID %s] \n%s" % (sequence.id,sequence.id,sequence.sequence)


if __name__=="__main__":
	main()
