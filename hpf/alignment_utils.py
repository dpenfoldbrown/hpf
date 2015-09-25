
from Bio.Align.Generic import Alignment

#kdrew: adapted from # http://biostumblematic.wordpress.com
def percent_identity(alignment):
	j = 0
	i = 0
	pi_dict = {}
	for record in alignment:
		for amino_acid in record.seq:
			if amino_acid == '-':
			    pass
			else:
			    if amino_acid == alignment[0].seq[j]:
				i += 1
			j += 1
		j = 0
		seq = str(record.seq)
		gap_strip = seq.replace('-', '')
		percent = 100*i/len(gap_strip)
		#print record.id+' '+str(percent)
		pi_dict[record.id] = percent
		i=0	

	print pi_dict
	return pi_dict
