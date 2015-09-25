
from hpf_report.model.meta import Session
from hpf.hddb.db import Domain, Protein, Sequence, DomainRegion, FilesystemOutfile, DomainFoldableMap
from itertools import imap
from sqlalchemy import not_
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm import aliased
import math
from operator import itemgetter


def main():

	#sequence_key = 271261
	#sequence_key = 93061
	#map_sequence_key(93061)
	#map_sequence_key(271261)
	for seq_key in Session().query(FilesystemOutfile.sequence_key).filter(FilesystemOutfile.sequence_key <> FilesystemOutfile.parent_sequence_key).distinct().all():
		print seq_key[0]
		map_sequence_key(seq_key[0])


def map_sequence_key(sequence_key):

	s1 = aliased(Sequence)
	s2 = aliased(Sequence)
	session = Session()
	sequences = session.query(FilesystemOutfile,s1,s2).join((s1,FilesystemOutfile.sequence), (s2,FilesystemOutfile.parent_sequence)).filter(FilesystemOutfile.sequence_key==sequence_key).all()

	for seq in sequences:

		print seq
		#print seq[0].parent_sequence_key
		fold_start = seq[2].sequence.find(seq[1].sequence) + 1
		fold_stop = fold_start + len(seq[1].sequence)
		print fold_start, fold_stop

		domain_regions = session.query(Domain, DomainRegion).join(DomainRegion.domain).filter(Domain.parent_sequence_key==seq[0].parent_sequence_key).all()
		for dreg in domain_regions:
			print dreg
			domain_start = dreg[1].start
			domain_stop = dreg[1].stop

			loc_dict = {"fold_start":fold_start, "fold_stop":fold_stop, "domain_start":domain_start, "domain_stop":domain_stop}
			coverage = percent_coverage(loc_dict)
			print coverage
			#kdrew: do not add entries which are nonoverlapping
			if not coverage['location'] == "fold_domain_nonoverlapping" and not coverage['location'] == "domain_fold_nonoverlapping":
				dfm = DomainFoldableMap(parent_sequence_key = seq[0].parent_sequence_key, fold_sequence_key = seq[0].sequence_key, domain_sequence_key = dreg[0].domain_sequence_key, domain_key = dreg[0].id, outfile_key = seq[0].id, fold_start=fold_start, fold_stop=fold_stop, domain_start=domain_start, domain_stop=domain_stop, fold_coverage=coverage['fold_coverage'], domain_coverage=coverage['domain_coverage'])
				print "dfm: ", dfm

				session.add_all([dfm,])

def percent_coverage(loc_dict):
	location = ""

	fold_len = loc_dict["fold_stop"] - loc_dict["fold_start"] + 1
	domain_len = loc_dict["domain_stop"] - loc_dict["domain_start"] + 1

	loc_sort = sorted(loc_dict.iteritems(), key=itemgetter(1))
	print loc_sort

	#kdrew: test to see if the second position is either end making them nonoverlapping
	if loc_sort[1][0] == 'fold_stop': 
		print "nonoverlapping: fold before domain"
		return {"fold_coverage":0,"domain_coverage":0,"location":"fold_domain_nonoverlapping"}
	if loc_sort[1][0] == 'domain_stop':
		print "nonoverlapping: fold after domain"
		return {"fold_coverage":0,"domain_coverage":0,"location":"domain_fold_nonoverlapping"}

	if loc_sort[1][0] == 'domain_start' and loc_sort[2][0] == 'fold_stop':
		location = "fold_span_domain_start"
	if loc_sort[1][0] == 'fold_start' and loc_sort[2][0] == 'fold_stop':
		location = "fold_encompassed"
	if loc_sort[1][0] == 'fold_start' and loc_sort[2][0] == 'domain_stop':
		location = "fold_span_domain_stop"
	if loc_sort[1][0] == 'domain_start' and loc_sort[2][0] == 'domain_stop':
		location = "domain_encompassed"

	#kdrew: subtract the values to get the amount of overlap
	num_overlap = loc_sort[2][1] - loc_sort[1][1] + 1
	return {"fold_coverage":num_overlap/(1.0*fold_len),"domain_coverage":num_overlap/(1.0*domain_len),"location":location}


if __name__ == "__main__":
	main()

