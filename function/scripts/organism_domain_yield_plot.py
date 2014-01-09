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
	print "Welcome to Organism Domain Yield"

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

	ody = OrganismDomainYield(config_file, config_section)
	organisms_dict, organisms_name  = ody.getYields()
	#ody.plotYieldsPie(name, yield_dict)
	ody.plotYieldsChart(organisms_dict,organisms_name, normalize=True)

	fr_sum = 0
	psiblast_sum = 0
	msa_sum = 0
	pfam_sum = 0
	heuristic_sum = 0
	for key in ody.experiment_keys:
		normalized_org = ody.normalize(organisms_dict[key])
		print organisms_name[key]
		for dkey in ody.domain_keys:
			print dkey, normalized_org[dkey]
		print "psiblast + fr: ", normalized_org["PDB-Blast"] + normalized_org["FFAS03"]
		fr_sum += normalized_org["FFAS03"]
		psiblast_sum += normalized_org["PDB-Blast"]
		msa_sum += normalized_org["MSA"]
		pfam_sum += normalized_org["Pfam"]
		heuristic_sum += normalized_org["Heuristic"]

	print "average fr: ", fr_sum/len(ody.experiment_keys)
	print "average psiblast: ", psiblast_sum/len(ody.experiment_keys)
	print "average msa: ", msa_sum/len(ody.experiment_keys)
	print "average pfam: ", pfam_sum/len(ody.experiment_keys)
	print "average heuristic: ", heuristic_sum/len(ody.experiment_keys)
		

class OrganismDomainYield():
	
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

		#self.experiment_key = config.get(section,'experiment_key')
		self.experiment_keys_str = config.get(section,'experiment_keys')
		self.experiment_keys = map(int,self.experiment_keys_str.split(','))
		print self.experiment_keys
		self.yield_plot_filename = config.get(section,'yield_plot_filename')
		
		self.domain_keys = ['PDB-Blast','FFAS03','MSA','Pfam','Heuristic']
		self.color = {'Pfam':'#757116','Heuristic':'#AEBC21','MSA':'#D9DB56','FFAS03':'#00477F','PDB-Blast':'#4C88BE'}
		self.ddb_conn = MySQLdb.connect(host=self.ddb_host, user=self.sql_user, passwd=self.sql_password, db=self.ddb_db, port=self.ddb_port)

	#kdrew: gets yield data
	def getYields(self):
		cursor = self.ddb_conn.cursor(MySQLdb.cursors.DictCursor)

		#kdrew: single organism
		#query = "select * from p_ginzu where id = %s"
		#cursor.execute(query,(self.experiment_key,))
		#row = cursor.fetchone ()

		#query = "select * from p_ginzu"
		#kdrew: old sample query result
		#+-----+------------+----------------------+---------------------+----------+------------------+------+-------+------------+
		#| id  | name       | ifnull(d.proteins,0) | ifnull(d.domains,0) | psiblast | fold_recognition | pfam | msa   | unassigned |
		#+-----+------------+----------------------+---------------------+----------+------------------+------+-------+------------+
		#| 804 | H. sapiens |                31801 |               79064 |    38047 |            10167 | 2028 | 12061 |      15365 | 
		#+-----+------------+----------------------+---------------------+----------+------------------+------+-------+------------+

		query = """SELECT pe.id, pe.name, pds.domain_type, SUM(pds.size) AS size 
				FROM p_domain_sccs AS pds, p_experiments AS pe 
				WHERE pe.id = pds.id 
				AND pe.id  = %s
				AND pds.domain_type IN ('pfam','fold_recognition', 'msa','unassigned','psiblast') 
				GROUP BY pe.id, pds.domain_type 
				ORDER BY pe.id, size;
			"""
		#kdrew: new sample query result
		#+-----+-----------------+------------------+-------+
		#| id  | name            | domain_type      | size  |
		#+-----+-----------------+------------------+-------+
		#| 804 | H. sapiens      | pfam             |  2013 | 
		#| 804 | H. sapiens      | fold_recognition | 10149 | 
		#| 804 | H. sapiens      | msa              | 12056 | 
		#| 804 | H. sapiens      | unassigned       | 15358 | 
		#| 804 | H. sapiens      | psiblast         | 37895 | 

		organism_yields = {}
		organism_name = {}
		for key in self.experiment_keys:
			cursor.execute(query, (key))
			rows = cursor.fetchall()
			organism_yields[key] = {}

			for row in rows:
				if row['domain_type'] == 'unassigned':
					organism_yields[key]['Heuristic'] = row['size']
				elif row['domain_type'] == 'psiblast':
					organism_yields[key]['PDB-Blast'] = row['size']
				elif row['domain_type'] == 'fold_recognition':
					organism_yields[key]['FFAS03'] = row['size']
				elif row['domain_type'] == 'msa':
					organism_yields[key]['MSA'] = row['size']
				elif row['domain_type'] == 'pfam':
					organism_yields[key]['Pfam'] = row['size']
				else:
					organism_yields[key][row['domain_type']] = row['size']
				organism_name[row['id']] = row['name']

		return organism_yields, organism_name

	def plotYieldsPie(self,name, yield_dict):
		pylab.figure(figsize=(10,10))
		ax = pylab.axes([0.3, 0.3, 0.4, 0.4])
		patches, texts = pylab.pie(yield_dict.values(),labels=yield_dict.keys())
		pylab.title("Ginzu domain frequencies: "+ name)
		#kdrew: a little filename manipulation
		f_parts = self.yield_plot_filename.rpartition('.')
		filename = f_parts[0] + "_" + self.experiment_key + "." + f_parts[2]
		pylab.savefig(filename)

	def plotYieldsChart(self,organisms_dict, organisms_name, normalize=True):
		#pylab.figure(figsize=(7,7))
		ax = pylab.axes([0.2, 0.1, .7, .7])
		ind = arange(len(organisms_dict))		
		width = .35
		plots = {}
		domain_dict = {'PDB-Blast':[],'FFAS03':[],'MSA':[],'Pfam':[],'Heuristic':[]}
		for key in self.experiment_keys:
			for dkey in self.domain_keys:
				try:
					if normalize:
						normalized_org = self.normalize(organisms_dict[key])
						domain_dict[dkey].append(normalized_org[dkey])
					else:
						domain_dict[dkey].append(organisms_dict[key][dkey])
				except KeyError:
					domain_dict[dkey].append(0)


		sum_bottom = [0 for i in arange(len(organisms_dict))]
		print "length of sum_bottom:", len(sum_bottom)
		for dkey in self.domain_keys:
			#plots[dkey] = pylab.bar(ind, domain_dict[dkey], color=self.color[dkey],bottom=sum_bottom)
			plots[dkey] = pylab.barh(ind, domain_dict[dkey], color=self.color[dkey],left=sum_bottom)
			sum_bottom = map(sum, zip(sum_bottom,domain_dict[dkey]))
			
		#pylab.xticks(ind+width/2, organisms_dict.keys(), rotation=45)
		pylab.yticks(ind+width/2, [organisms_name[i] for i in self.experiment_keys])
		pylab.xticks()
		pylab.title("Domain type coverage")
		
		fp = font.FontProperties(size="x-small")
		pylab.legend( [plots[key][0] for key in self.domain_keys], self.domain_keys, prop=fp, markerscale=.5, loc=(.85,.85))


		#kdrew: a little filename manipulation
		#f_parts = self.yield_plot_filename.rpartition('.')
		#filename = f_parts[0] + "_" + self.experiment_key + "." + f_parts[2]
		pylab.savefig(self.yield_plot_filename)

	def normalize(self,yield_dict):
		y_dict = {}
		for key in yield_dict.keys():
			y_dict[key] = yield_dict[key]/sum(yield_dict.values())

		return y_dict

if __name__ == "__main__":
	main()
    


