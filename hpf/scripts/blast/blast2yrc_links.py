
import unittest
import os
import sys
import MySQLdb
import getopt
import ConfigParser
import datetime

from Bio.Blast import NCBIXML
from hpf.function.createdbs.blast_filter import Filtered, BlastFilter

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to blast2yrc_links"

	config_file = "config.cfg"
	config_section = "Default"

	try:
		print "before getopt"
		opts, args = getopt.getopt(sys.argv[1:], "hc:s:", ["help", "config=", "section="]) 
		print "after getopt"
		print opts
		print args
	except getopt.GetoptError:	
		usage()
		sys.exit(2)
	for opt, arg in opts:
		print opt
		print arg
		if opt in ("-h","--help"):
			usage()
			sys.exit()
		elif opt in ("-c","--config"):
			config_file = arg
		elif opt in ("-s","--section"):
			config_section = arg
		else:
			assert False, "unhandled option"

	bsk = Blast2YRC(config_file, config_section)
	filtered = bsk.runblast2seq_key()
	print "filtered length: ",len(filtered)

	#kdrew: from pwinters 2010-11-09 email
	# the database id's are our sequence keys
	# YRC maps to the table hpf.yeastrcu ... this is the sequence -> sequence -> proteinID mapping
	# YRCProtein maps to the protein names and descriptors for each yrc proteinID
	from hpf.hddb.db import YRC, YRCProtein, Domain, DomainRegion, DomainSCCS, Session
	session = Session()
	for filt in filtered:
		yrc = session.query(YRC).filter(YRC.sequence_key==int(filtered[filt].hit_id)).all()
		links = ["http://www.yeastrc.org/pdr/viewProtein.do?id=%i" % y.yrc_protein_key for y in yrc]
		#print links
		link = ""

		for y in yrc:
			#print y
			yrc_protein_entry = session.query(YRCProtein).filter(YRCProtein.yrc_protein_key==y.yrc_protein_key).all()
			#kdrew: this is kinda lame and does not work for all cases, currently tests to make sure the taxid is not 0 but if yrc_protein_key not in YRCProtein table just skips test
			#kdrew: need to update YRCProtein table so that it includes all protein ids
			#kdrew: might be useful to have a list of taxids to compare against
			if len(yrc_protein_entry) > 0:
				#print yrc_protein_entry[0]
				if yrc_protein_entry[0].taxonomy_id == 0:
					#kdrew: no taxid so move to the next yrc_protein_key
					continue
				else:
					link = "http://www.yeastrc.org/pdr/viewProtein.do?id=%i" % y.yrc_protein_key
					#print link
					#kdrew: got a valid yrc_protein_key with taxid, break out of loop
					break
			else:
				link = "http://www.yeastrc.org/pdr/viewProtein.do?id=%i" % y.yrc_protein_key
				#print link
				#kdrew: yrc_protein_key not in YRCProtein so use it regardless, break out of loop
				break

		print filtered[filt].query_id, "|", filtered[filt].hit_id, "|", link

		domains = session.query(DomainSCCS).filter(DomainSCCS.parent_sequence_key == int(filtered[filt].hit_id)).all()
		for domain in domains:
			domain_entry = session.query(Domain).filter(Domain.domain_sequence_key == domain.domain_sequence_key).first()
			domainregion = session.query(DomainRegion).filter(DomainRegion.domain_key == domain_entry.id).first()
			
			print "|","|","|",domain.domain_type,"|",domain.sccs,"|",domain.pdbid,"|",domain.chain,"|",domainregion.start,"|",domainregion.stop,"|",domain.confidence

class Blast2YRC():
	
	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		self.e_value_threshold = config.getfloat(section, 'e_value_threshold')
		self.length_threshold = config.getfloat(section, 'length_threshold')
		self.identity_cutoff = config.getfloat(section, 'identity_cutoff')

		self.blast_file = config.get(section, 'blast_file')

	def runblast2seq_key(self):
		outfile_handle = open(self.blast_file)
		blast_records = NCBIXML.parse(outfile_handle)
		ba = PSSPBlastFilter(eval_cutoff = self.e_value_threshold, length_cutoff = self.length_threshold, identity_cutoff = self.identity_cutoff)

		filtered = ba.filterBlast(blast_records)

		print len(filtered)
		for filt in filtered:
			print filtered[filt].query_id, filtered[filt].hit_id

		return filtered

class QualityFiltered(Filtered):
	def __init__(self,q_id,h_id,h_sf):
		Filtered.__init__(self,q_id,h_id)
		self.superfamily = h_sf

class PSSPBlastFilter(BlastFilter):
	#kdrew: overriden to handle specific parsing
	def parse_blast_hit(self, record_query, align_def,hit_title=None):
		print "rq:", record_query
		print "ad:", align_def
		print "ht:", hit_title
		query_id = record_query.split('|')[1].strip()
		hit_id = hit_title.split('|')[1].split()[0].strip()

		return Filtered(query_id, hit_id)



if __name__ == "__main__":
	main()


