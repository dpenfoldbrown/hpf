from __future__ import division

import unittest
import os
import sys
import MySQLdb
import getopt
import ConfigParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import matplotlib.font_manager as font

import networkx as nx
#kdrew: this is so it creates the figures without accessing the windowing system

from itertools import izip

import pylab
from numpy import arange

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Domain Count"

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

	dc = DomainCount(config_file, config_section)
	domain_size, domain_counts  = dc.getCount()
	dc.plotDomainChart(domain_size,domain_counts)
	dc.getExpNames()

	weight_sum = 0
	j = 0
	for i in domain_size:
		weight_sum += i*domain_counts[j]
		j+=1

	print "Average domain size: ", weight_sum/sum(domain_counts)

class DomainCount():
	
	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		#kdrew: blast parameters
		self.ddb_db = config.get(section, 'ddb_db')
		self.ddb_host= config.get(section, 'ddb_host')
		self.ddb_port = config.getint(section, 'ddb_port')
		print self.ddb_port

		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')

		self.experiment_keys_str = config.get(section,'experiment_keys')
		if self.experiment_keys_str != "":
			self.experiment_keys = map(int,self.experiment_keys_str.split(','))
		print self.experiment_keys_str

		self.plot_filename = config.get(section,'plot_filename')
		
		#self.domain_keys = ['psiblast','fold_recognition','msa','pfam','heuristic']
		#self.color = {'pfam':'#757116','heuristic':'#AEBC21','msa':'#D9DB56','fold_recognition':'#00477F','psiblast':'#4C88BE'}

		self.ddb_conn = MySQLdb.connect(host=self.ddb_host, user=self.sql_user, passwd=self.sql_password, db=self.ddb_db, port=self.ddb_port)

	def getExpNames(self):
		cursor = self.ddb_conn.cursor(MySQLdb.cursors.DictCursor)
		query ="""
				SELECT name 
				FROM p_experiments AS pe
				WHERE pe.id in ("""+self.experiment_keys_str+""")
			"""
		cursor.execute(query)
		rows = cursor.fetchall()

		for row in rows:
			print row['name']
		

	#kdrew: gets yield data
	def getCount(self):
		cursor = self.ddb_conn.cursor(MySQLdb.cursors.DictCursor)

		exp_keys = ""
		if self.experiment_keys_str != "":
			exp_keys = " WHERE pds.id in ("+self.experiment_keys_str+")"

		query = """
				SELECT size, COUNT(distinct parent_sequence_key) AS cnt 
				FROM (SELECT parent_sequence_key, SUM(size) AS size 
					FROM p_domain_sccs as pds
			"""+exp_keys+"""
					GROUP BY parent_sequence_key 
					ORDER BY size DESC) 
					AS tmp 
				GROUP BY size 
				ORDER BY cnt DESC
			"""

		domain_size= []
		domain_count = [] 

		cursor.execute(query) 
		rows = cursor.fetchall()

		for row in rows:
			domain_size.append(int(row['size']))
			domain_count.append(int(row['cnt']))

		print domain_size, domain_count
		return domain_size, domain_count

	def plotDomainChart(self, domain_size, domain_count):
		pylab.bar(domain_size, domain_count)

		pylab.xticks()
		pylab.title("Domains per protein")
		pylab.xlabel("Num. of domains")
		pylab.ylabel("Frequency")
		
		pylab.savefig(self.plot_filename)


if __name__ == "__main__":
	main()
    


