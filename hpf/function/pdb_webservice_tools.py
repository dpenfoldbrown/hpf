
from cStringIO  import StringIO
from datetime import date
import SOAPpy
from Bio.Blast import NCBIXML

from hpf.function.createdbs.blast_filter import Filtered, BlastFilter

class PDBBlastFilter(BlastFilter):
	#kdrew: overriden to handle specific parsing
	def parse_blast_hit(self, record_query, align_def):
		pdb_id = align_def.split(":")[0]
		#print pdb_id
		return Filtered(record_query, pdb_id)

PDB_WS_URL = "http://www.pdb.org/pdb/services/pdbws"
server = SOAPpy.SOAPProxy(PDB_WS_URL)

#kdrew: this function uses pdb's webservice to blast, then filters by specified thresholds and then uses the webservice to get release date
def get_pdb_date_Blast(sequence, e_value_threshold, length_threshold, identity_cutoff):

	blastOutput = server.blastPDB ( sequence, e_value_threshold, "BLOSUM62", "XML" );
	blast_out = StringIO(blastOutput)
	#print blastOutput

	blast_records = NCBIXML.parse(blast_out)
	ba = PDBBlastFilter(eval_cutoff = e_value_threshold, length_cutoff = length_threshold, identity_cutoff = identity_cutoff, multi_hits=True)

	filtered = ba.filterBlast(blast_records)

	return get_pdb_date(filtered)


def get_pdb_date(filtered_one,parse_hit_id = False):
	date_dict = {}
	for hit in filtered_one:
		h_id = hit.hit_id
		if parse_hit_id:
			h_id = hit.hit_id[0:4]
			#print h_id
		releaseDates = server.getReleaseDates( (h_id,))
		cdate = convertDate(releaseDates[0].split()[1])
		date_dict[cdate] = releaseDates[0].split()[0]

	return date_dict
def get_min_date(date_dict):
	if len(date_dict) > 0:
		return min(date_dict) , date_dict[min(date_dict)]
	else:
		print "no dates in dictionary"
		return None

def convertDate(iso_string):
	return date(int(iso_string.split("-")[0]),int(iso_string.split("-")[1]),int(iso_string.split("-")[2]))
	 
def main():
	e_value_threshold = 1e-2
	length_threshold = .5
	identity_cutoff = 0
	sequence = "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNAL"
	pdb_dates = get_pdb_date_Blast(sequence, e_value_threshold, length_threshold, identity_cutoff)
	print get_min_date(pdb_dates)

if __name__ == "__main__":
	main()

