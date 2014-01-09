
import unittest
import os
import MySQLdb
import getopt
import ConfigParser
import shelve
import matplotlib
#kdrew: this is so it creates the figures without accessing the windowing system
matplotlib.use('Agg')
import matplotlib.pyplot
import matplotlib.mlab as mlab
import networkx as nx
from latex_format import LatexFormatter


import pylab
from numpy import arange
import numpy

from hpf.function.metric import *
from hpf.function.structure import *
from hpf.function.domain import *
from hpf.function.bayes import *
from hpf.function.term import *
from hpf.function.createdbs.astral import *
from hpf.go import GODAG,_pydot

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Since Solved Analysis"

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

	ds = DenovoStats(config_file, config_section)

	ds.runDenovoStats()


class DenovoStats():
	
	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		#kdrew: blast parameters
		self.hpf_db = config.get(section, 'hpf_db')
		self.hpf_host = config.get(section, 'hpf_host')
		self.hpf_port = config.getint(section, 'hpf_port')

		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')

		self.denovo_freq_plot_filename = config.get(section,'denovo_freq_plot_filename')
		self.high_conf_limit = config.getfloat(section,'high_conf_limit')

		self.by_random = config.getboolean(section, 'by_random')
		self.rand_seed = config.getint(section, 'rand_seed')
		self.sample_size = config.getint(section, 'sample_size')

		self.psiblast_flag = config.getboolean(section, 'psiblast_flag')
		self.fr_flag = config.getboolean(section, 'fr_flag')

		self.classA_color = "blue"
		self.classB_color = "orange"
		self.classC_color = "red"
		self.classD_color = "green"

		self.hpf_conn = MySQLdb.connect(host=self.hpf_host, user=self.sql_user, passwd=self.sql_password, db=self.hpf_db, port=self.hpf_port)

	def runDenovoStats(self):
		cursor = self.hpf_conn.cursor(MySQLdb.cursors.DictCursor)

		rand_sub_query = ""
		if self.by_random:
			tmp_create_query = """CREATE TEMPORARY TABLE tmp_rand SELECT DISTINCT domain_sequence_key as sequence_key FROM p_domain_sccs ORDER BY RAND(%s) LIMIT %s """
			cursor.execute(tmp_create_query, (self.rand_seed, self.sample_size))
			cursor.execute("ALTER TABLE tmp_rand ADD INDEX(sequence_key)")
			rand_sub_query = """ AND pds.domain_sequence_key IN (SELECT * FROM tmp_rand) """

		if self.psiblast_flag:
			where_sub_query = "WHERE pds.domain_type IN ('psiblast')"
		elif self.fr_flag:
			where_sub_query = "WHERE pds.domain_type IN ('fold_recognition')"
		else:
			where_sub_query = "WHERE pds.domain_type IN ('unassigned','pfam','msa')"


		#kdrew: new query from pwinters 12/1/09
		query = """
				SELECT SUBSTRING_INDEX(pds.sccs, '.',1) AS class, SUM(size) AS cnt 
				FROM p_domain_sccs AS pds 
				""" + where_sub_query + """
				AND pds.sccs IS NOT NULL AND pds.confidence > %s
				""" + rand_sub_query + """
				GROUP BY SUBSTRING_INDEX(pds.sccs,'.',1) ORDER BY SUBSTRING_INDEX(pds.sccs,'.',1)
			"""

		print query
		cursor.execute(query,(0,))
		all_rows = cursor.fetchall ()
	
		prob_true_list = list()
		prob_false_list = list()
		prob_list = list()

		t_dict = {}
		import string
		ltrs = string.lowercase
		#kdrew: initialize letters a-i to be zero
		for i in range(0,9):
			t_dict[ltrs[i]+"_total"] = 0
			t_dict[ltrs[i]+"_high_conf"] = 0
		
		for all_row in all_rows:
			print all_row["class"], all_row["cnt"]
			t_dict[all_row["class"]+"_total"] = all_row["cnt"]
			

		cursor.execute(query,(self.high_conf_limit,))
		high_conf_rows = cursor.fetchall ()

		for high_conf_row in high_conf_rows:
			print high_conf_row["class"], high_conf_row["cnt"]
			t_dict[high_conf_row["class"]+"_high_conf"] = high_conf_row["cnt"]
			

		lf = LatexFormatter()
		print "De novo byClass Latex table format:\n"
		print lf.denovo_byClass_table_format(t_dict)

	def plot_SS_Freq_Analysis(self,list1, c='b', label="", linetype='-', bins=10):
		n, bins, patches = pylab.hist(list1,visible=False,bins=bins)
		print n
		print bins
		pylab.xlim(xmin=.3, xmax=1.0)
		l = pylab.plot(bins[1:], n, linetype, color=c,label=label,linewidth=1)
		return n, bins



if __name__ == "__main__":
	main()


