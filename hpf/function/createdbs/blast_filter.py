
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML
from Bio import SeqIO


#kdrew: holds filtered blast records
class Filtered():
	def __init__(self,q_id, h_id):
		self.query_id = q_id
		self.hit_id = h_id

class BlastFilter():
	def __init__(self, blast_exe=None, blast_db=None, blast_query=None, eval_cutoff = 1e-8, length_cutoff = .85, identity_cutoff = 0, blast_prog = "blastp", blast_processors=1, multi_hits=False):
		#kdrew: location of executable
		self.blast_exe = blast_exe
		self.blast_db = blast_db
		#kdrew: file to be blasted
		self.blast_query = blast_query
		#kdrew: type of blast to be run (ie. blastp)
		self.blast_prog = blast_prog
		self.eval_cutoff = eval_cutoff
		self.length_cutoff = length_cutoff
		self.identity_cutoff = identity_cutoff
		self.blast_processors = blast_processors
		#kdrew: return just the top hit that passes the filter or all hits that pass filter
		self.return_multiple_hits = multi_hits

	#kdrew: run blast and threshold, return go seq_id with matched superfamily
	def runBlast(self, result_handle=None):
		if result_handle == None:
			result_handle, error_handle = NCBIStandalone.blastall(self.blast_exe, self.blast_prog, self.blast_db, self.blast_query, nprocessors=self.blast_processors)

		#kdrew: if we want to pre-run blast, just run this line
		blast_records = NCBIXML.parse(result_handle)
		return blast_records


	#kdrew: returns a list of blast matches which pass the filter
	def filterBlast(self, blast_records):
		filter_dict= {}
		for blast_record in blast_records:
			filtered_record = self.filterBlastOne(blast_record)
			if None != filtered_record:
				if self.return_multiple_hits:
					filter_dict[str(filtered_record[0].query_id)] = filtered_record
				else:
					filter_dict[str(filtered_record.query_id)] = filtered_record

		return filter_dict

	#kdrew: default but should be overridden
	def parse_blast_hit(self, record_query, align_def):
		query_id = record_query.split('|')[0].strip()
		hit_id = align_def.split()[0].strip()

		return Filtered(query_id, hit_id)

	def filterBlastOne(self, blast_record):
		filter_rec = None
		filter_rec_list = list()
		recorded_flag = False
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < self.eval_cutoff and hsp.align_length > alignment.length*self.length_cutoff and hsp.identities > self.identity_cutoff:
					#print "passed filter"
					filter_rec = self.parse_blast_hit(blast_record.query, alignment.hit_def)
					recorded_flag = True
					filter_rec_list.append(filter_rec)

			if not self.return_multiple_hits:
				#kdrew: if this is set it means we have already recorded a hit for this sequence
				if recorded_flag:
					break

		if self.return_multiple_hits:
			if 0 < len(filter_rec_list):
				return filter_rec_list
			else:
				return None
		else:
			return filter_rec


	def writeFasta(self, filter_dict, out_file=None):
		#print self.blast_query
		#print filter_dict

		handle_in = open(self.blast_query)
		if None == out_file:
			handle_out = open(self.blast_query.rpartition('.')[0]+"_filtered.fasta", "w")
		else:
			handle_out = open(out_file, "w")
		out_list = []
		for seq_record in SeqIO.parse(handle_in, "fasta") :
			#print seq_record.id
			#print repr(seq_record.seq)
			#print len(seq_record)
			if seq_record.id.split('|')[0] in filter_dict:
				print seq_record.id
				out_list.append(seq_record)

		SeqIO.write(out_list, handle_out, "fasta")
			
		handle_in.close()
		handle_out.close()


class AstralFiltered(Filtered):
	def __init__(self,q_id,h_id,h_sf):
		Filtered.__init__(self,q_id,h_id)
		self.superfamily = h_sf

class AstralBlastFilter(BlastFilter):
	#kdrew: overriden to handle specific parsing
	def parse_blast_hit(self, record_query, align_def):
		query_id = record_query.split('|')[0].strip()
		hit_id = align_def.split()[0].strip()
		sccs_id = align_def.split()[1].strip()
		superfamily = sccs_id.rpartition('.')[0]

		return AstralFiltered(query_id, hit_id, superfamily)

class GOBlastFilter(BlastFilter):
	#kdrew: overriden to handle specific parsing
	def parse_blast_hit(self, record_query, align_def):
		print "query_id_raw: ", record_query
		query_id = record_query.split('|')[1].strip()
		print "hit_id_raw: ", align_def
		hit_id = align_def.split('|')[0].strip()
		print query_id, hit_id

		return Filtered(int(query_id), int(hit_id))


