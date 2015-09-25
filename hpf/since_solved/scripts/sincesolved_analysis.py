
from __future__ import division
import operator

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
from hpf.function.scripts.latex_format import LatexFormatter


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

	ssa = SinceSolved_Analysis(config_file, config_section)
	if ssa.do_sincesolved_mammoth_analysis:
		ssa.runSSMammothAnalysis()
	if ssa.do_sincesolved_freq_analysis:
		if ssa.by_class_flag:
			if ssa.cor_VS_incor_flag:
				a_n,a_n_corr,b_n,b_n_corr,c_n,c_n_corr,d_n,d_n_corr,o_n,o_n_corr = ssa.plot_cor_VS_incor_by_class()
			else:
				a_n,a_n_corr,b_n,b_n_corr,c_n,c_n_corr,d_n,d_n_corr,o_n,o_n_corr = ssa.plot_by_class()
				lf = LatexFormatter()
				print "De novo byClass Latex table format:\n"
				print lf.denovo_byClass_table(a_n,a_n_corr,b_n,b_n_corr,c_n,c_n_corr,d_n,d_n_corr,o_n,o_n_corr)
		else:
			ssa.runSS_functionCorrelation_Analysis()

	if ssa.do_sincesolved_correlation_analysis or ssa.do_sincesolved_zscore_corr_analysis:
		if ssa.by_class_flag:
			ssa.by_class = 'a'
			plist, clist, zlist = ssa.runSS_functionCorrelation_Analysis()
			ssa.plot_SS_functionCorrelation_Analysis(zlist,clist,c=ssa.classA_color, label="SCOP Class A")
			zarray = numpy.array(zlist)
			carray = numpy.array(clist)
			m,c = ssa.lsqfit(zarray,carray)
			pylab.plot(zarray, m*zarray+c,'r',c=ssa.classA_color)

			ssa.by_class = 'b'
			plist, clist, zlist = ssa.runSS_functionCorrelation_Analysis()
			ssa.plot_SS_functionCorrelation_Analysis(zlist,clist,c=ssa.classB_color, label="SCOP Class B")
			zarray = numpy.array(zlist)
			carray = numpy.array(clist)
			m,c = ssa.lsqfit(zarray,carray)
			pylab.plot(zarray, m*zarray+c,'r',c=ssa.classB_color)

			ssa.by_class = 'c'
			plist, clist, zlist = ssa.runSS_functionCorrelation_Analysis()
			ssa.plot_SS_functionCorrelation_Analysis(zlist,clist,c=ssa.classC_color, label="SCOP Class C")
			zarray = numpy.array(zlist)
			carray = numpy.array(clist)
			m,c = ssa.lsqfit(zarray,carray)
			pylab.plot(zarray, m*zarray+c,'r',c=ssa.classC_color)


			ssa.by_class = 'd'
			plist, clist, zlist = ssa.runSS_functionCorrelation_Analysis()
			ssa.plot_SS_functionCorrelation_Analysis(zlist,clist,c=ssa.classD_color, label="SCOP Class D")
			zarray = numpy.array(zlist)
			carray = numpy.array(clist)
			m,c = ssa.lsqfit(zarray,carray)
			pylab.plot(zarray, m*zarray+c,'r',c=ssa.classD_color)


			pylab.suptitle("Incorrect predicted superfamilies vs true superfamilies (mcm score > "+str(ssa.mcm_value_limit)+")")
			pylab.xlabel("Mammoth Z-score\n(measure of structural similarity)")
			pylab.ylabel("GO function correlation")
			leg = pylab.legend()
			for t in leg.get_texts():
				t.set_fontsize('small')
			pylab.savefig(ssa.sincesolved_zscore_corr_plot_filename)
		else:
			ssa.runSS_functionCorrelation_Analysis()
	

class SinceSolved_Analysis():
	
	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		#kdrew: blast parameters
		self.metric_db = config.get(section, 'metric_db')
		self.metric_host= config.get(section, 'metric_host')
		self.mygo_db = config.get(section, 'mygo_db')
		self.mygo_host= config.get(section, 'mygo_host')
		self.since_solved_db = config.get(section, 'since_solved_db')
		self.since_solved_host = config.get(section, 'since_solved_host')

		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')
		self.prob_table_name = config.get(section, 'probability_table_name')

		self.noLoad = config.getboolean(section,'noLoad')
		self.probability_shelve = config.get(section,'probability_shelve')

		self.do_sincesolved_mammoth_analysis = config.getboolean(section,'do_sincesolved_mammoth_analysis')
		self.sincesolved_mammoth_plot_filename = config.get(section,'sincesolved_mammoth_plot_filename')

		self.do_sincesolved_correlation_analysis = config.getboolean(section,'do_sincesolved_correlation_analysis')
		self.sincesolved_correlation_plot_filename = config.get(section,'sincesolved_correlation_plot_filename')

		self.do_sincesolved_zscore_corr_analysis = config.getboolean(section,'do_sincesolved_zscore_corr_analysis')
		self.sincesolved_zscore_corr_plot_filename = config.get(section,'sincesolved_zscore_corr_plot_filename')

		self.do_sincesolved_freq_analysis= config.getboolean(section,'do_sincesolved_freq_analysis')
		self.sincesolved_freq_plot_filename = config.get(section,'sincesolved_freq_plot_filename')

		self.intersection_limit = config.getint(section,'intersection_limit')
		self.mcm_value_limit = config.getfloat(section,'mcm_value_limit')
		self.scop_level = config.getint(section,'scop_level')
		self.plot_bins_str = config.get(section,'plot_bins')
		if self.plot_bins_str == "":
			self.plot_bins = 10
		else:
			self.plot_bins = map(float,list(self.plot_bins_str.split(',')))
		self.by_class_flag = config.getboolean(section,'by_class_flag')

		self.cor_VS_incor_flag = config.getboolean(section,'cor_VS_incor_flag')
		self.by_true_class_flag = config.getboolean(section,'by_true_class_flag')
		self.group_sf_flag = config.getboolean(section,'group_sf_flag')
		self.by_class = ""
		self.date_cutoff = config.get(section,'date_cutoff')
		self.redux = config.getboolean(section,'redux')
		self.paper = config.getboolean(section,'paper')
		self.over_flag = config.getboolean(section,'over_flag')
		self.under_flag = config.getboolean(section,'under_flag')
		self.over_under_value = config.getint(section,'over_under_value')

		self.exclude_clusters_str = config.get(section,'exclude_clusters')
		self.exclude_clusters = []
		if self.exclude_clusters_str != "":
			self.exclude_clusters = map(int, list(self.exclude_clusters_str.split(',')))


		self.classA_color = "blue"
		self.classB_color = "orange"
		self.classC_color = "red"
		self.classD_color = "green"

		dict = shelve.open(self.probability_shelve)
		
		self.prob_metric = Probability(dict=dict)

		self.mygo_conn = MySQLdb.connect(host=self.mygo_host, user=self.sql_user, passwd=self.sql_password, db=self.mygo_db)

		if not self.noLoad:
			self.metric_conn = MySQLdb.connect(host=self.metric_host, user=self.sql_user, passwd=self.sql_password, db=self.metric_db)
			self.prob_metric.load_metric(self.metric_conn, self.prob_table_name)

		self.since_solved_conn = MySQLdb.connect(host=self.since_solved_host, user=self.sql_user, passwd=self.sql_password, db=self.since_solved_db)

		if self.over_flag or self.under_flag:
			self.exclude_overunder_clusters()

	def exclude_overunder_clusters(self):
		cursor = self.since_solved_conn.cursor(MySQLdb.cursors.DictCursor)
		if self.over_flag:
			query = "SELECT distinct cluster_id from since_solved.clstr_seq_lenavg where avg_seq_len < %s"
		elif self.under_flag:
			query = "SELECT distinct cluster_id from since_solved.clstr_seq_lenavg where avg_seq_len >= %s"

		cursor.execute(query,(self.over_under_value,))
		rows = cursor.fetchall ()
		for row in rows:
			self.exclude_clusters.append(row['cluster_id'])

	def runSSFreqAnalysis(self):
		cursor = self.since_solved_conn.cursor(MySQLdb.cursors.DictCursor)
		if self.by_class != "":
			if self.by_true_class_flag:
				by_class_query = " and SUBSTRING_INDEX(ss.sccs,'.',1) = '"+self.by_class+"'"
				by_class_query2 = by_class_query
				if self.paper:
					by_class_query = " and SUBSTRING_INDEX(true_sf,'.',1) = '"+self.by_class+"'"
			else:
				by_class_query = " and SUBSTRING_INDEX(mc.experiment_sccs,'.',1) = '"+self.by_class+"'"
				by_class_query2 = " and SUBSTRING_INDEX(mcr.experiment_sccs,'.',1) = '"+self.by_class+"'"
				if self.paper:
					by_class_query = " and SUBSTRING_INDEX(pred_sf,'.',1) = '"+self.by_class+"'"
					if self.by_class == 'o':
						by_class_query = " and SUBSTRING_INDEX(pred_sf,'.',1) not in ('a','b','c','d')"
		else:
			by_class_query = ""
			by_class_query2 = ""

		if self.group_sf_flag:
			group_sf_query = " GROUP BY SUBSTRING_INDEX(mc.experiment_sccs,'.',%s) " % (self.scop_level,)
			group_sf_query2 = " GROUP BY SUBSTRING_INDEX(mcr.experiment_sccs,'.',%s) " % (self.scop_level,)
		else:
			group_sf_query = ""
			group_sf_query2 = ""

		if self.paper:
					#SELECT * from since_solved.p_mcm_cluster WHERE substring_index(true_sf,'.',1) in ('a','b','c','d')
			query = """ 
					SELECT * from since_solved.p_mcm_cluster WHERE 1=1
					""" + by_class_query + """ 
					ORDER BY cluster_id 
				"""
		else:
			query = """
					SELECT mc.cluster_id AS cluster_id, SUBSTRING_INDEX(mc.experiment_sccs,'.',%s) AS pred_sf, SUBSTRING_INDEX(ss.sccs,'.',%s) AS true_sf, mc.probability as top_prob
					FROM since_solved.mcm_cluster mc, since_solved.since_solved ss,
						(SELECT cluster_id, MAX(probability) as max_prob from since_solved.mcm_cluster group by cluster_id) as max_table 
					WHERE mc.cluster_id = ss.cluster_id and max_table.cluster_id = mc.cluster_id and mc.probability = max_table.max_prob 
						AND ss.ascession_date >= DATE(%s)
						AND mc.probability > %s
				""" + by_class_query + group_sf_query + """ 
					ORDER BY cluster_id
				"""

		if self.redux:
			if self.paper:
						#SELECT * from since_solved.p_mcm_redux WHERE substring_index(true_sf,'.',1) in ('a','b','c','d')
				query = """ 
						SELECT * from since_solved.p_mcm_redux WHERE 1=1
						""" + by_class_query +""" 
						ORDER BY cluster_id 
					"""
					
			else:
				query = """
						SELECT distinct cc.cluster_id AS cluster_id, SUBSTRING_INDEX(mcr.experiment_sccs,'.',%s) AS pred_sf, SUBSTRING_INDEX(ss.sccs,'.',%s) AS true_sf, mcr.probability as top_prob
						FROM since_solved.cdhit_clstr cc, since_solved.since_solved ss, since_solved.mcm_redux mcr,
							(SELECT cc.cluster_id, MAX(mcr.probability) AS max_prob FROM since_solved.cdhit_clstr cc, since_solved.mcm_redux as mcr WHERE mcr.sequence_key = cc.foreign_key group by cc.cluster_id) as max_table 
						WHERE cc.cluster_id = ss.cluster_id and max_table.cluster_id = cc.cluster_id and mcr.probability = max_table.max_prob and mcr.sequence_key = cc.foreign_key
							AND ss.ascession_date >= DATE(%s)
							AND mcr.probability > %s
					""" + by_class_query2 + group_sf_query2 + """ 
						ORDER BY cluster_id 
					"""


		print query
		if self.paper:
			cursor.execute(query,)
		else:
			cursor.execute(query,(self.scop_level,self.scop_level,self.date_cutoff,self.mcm_value_limit,))
		rows = cursor.fetchall ()
	
		prob_true_list = list()
		prob_false_list = list()
		prob_list = list()
		
		for row in rows:
			if int(row["cluster_id"]) not in self.exclude_clusters:
				print row["cluster_id"], row["pred_sf"], row["true_sf"], row["top_prob"]
				if row["pred_sf"] == row["true_sf"]:
					prob_true_list.append(row["top_prob"])
				else:
					prob_false_list.append(row["top_prob"])

				prob_list.append(row["top_prob"])
			#else:
			#	print "excluded:", row["cluster_id"], row["pred_sf"], row["true_sf"], row["top_prob"]
			
			
		print "found: ", len(prob_list)
		return prob_list, prob_true_list, prob_false_list

	def plot_SS_Freq_Analysis(self,list1, c='b', label="", linetype='-', marker='', markersize=4, bins=10,visible=True,ymax=125):
		print len(list1)
		n, bins, patches = pylab.hist(list1,visible=False,bins=bins)
		print n
		print bins
		pylab.xlim(xmin=.3, xmax=1.005)
		pylab.ylim(ymax=ymax)
		l = pylab.plot(bins[1:], n, linetype, marker=marker, markersize=markersize, color=c,label=label,linewidth=1, visible=visible)
		return n, bins


	def runSSMammothAnalysis(self):
		cursor = self.since_solved_conn.cursor(MySQLdb.cursors.DictCursor)

		query = """
				SELECT mc.cluster_id AS cluster_id, SUBSTRING_INDEX(mc.experiment_sccs,'.',3) AS pred_sf, SUBSTRING_INDEX(ss.sccs,'.',3) AS true_sf, AVG(mc.probability) AS avg_prob, m.* 
				FROM since_solved.mcm_cluster mc, since_solved.since_solved ss, pdb.mammoth as m 
				WHERE ((mc.cluster_id = ss.cluster_id) and m.prediction_sccs = SUBSTRING_INDEX(mc.experiment_sccs,'.',3) and SUBSTRING_INDEX(ss.sccs,'.',3) = m.experiment_sccs) 
					AND ss.ascession_date >= DATE(%s)
				GROUP BY mc.cluster_id, SUBSTRING_INDEX(mc.experiment_sccs,'.',3)
			"""

		cursor.execute(query, (self.date_cutoff,))
		rows = cursor.fetchall ()
	
		prob_true_list = list()
		zscore_true_list = list()
		prob_false_list = list()
		zscore_false_list = list()
		
		for row in rows:
			if row["pred_sf"] == row["true_sf"]:
				prob_true_list.append(row["avg_prob"])
				zscore_true_list.append(row["zscore"])
			else:
				prob_false_list.append(row["avg_prob"])
				zscore_false_list.append(row["zscore"])
			
		pylab.scatter(prob_false_list, zscore_false_list)
		pylab.savefig(self.sincesolved_mammoth_plot_filename)

	def runSS_functionCorrelation_Analysis(self):
		cursor = self.since_solved_conn.cursor(MySQLdb.cursors.DictCursor)
		
		if self.by_class != "":
			by_class_query = " and SUBSTRING_INDEX(mc.experiment_sccs,'.',1) = '"+self.by_class+"'"
		else:
			by_class_query = ""

		query = """
				SELECT mc.cluster_id AS cluster_id,SUBSTRING_INDEX(mc.experiment_sccs,'.',3) AS pred_sf, SUBSTRING_INDEX(ss.sccs,'.',3) AS true_sf, AVG(mc.probability) AS avg_prob, c.*, m.*
				FROM since_solved.mcm_cluster mc, since_solved.since_solved ss, pdb.correlation as c, pdb.mammoth as m
				WHERE ((mc.cluster_id = ss.cluster_id) and c.sf1 = SUBSTRING_INDEX(mc.experiment_sccs,'.',3) 
					and SUBSTRING_INDEX(ss.sccs,'.',3) = c.sf2 and c.sf1 = m.prediction_sccs and c.sf2 = m.experiment_sccs) 
					AND ss.ascession_date >= DATE(%s)
					and mc.probability > %s
			""" + by_class_query + """
				GROUP BY mc.cluster_id, SUBSTRING_INDEX(mc.experiment_sccs,'.',3)
			"""

		print query
		cursor.execute(query, (self.date_cutoff, self.mcm_value_limit,))
		rows = cursor.fetchall ()
	
		prob_true_list = list()
		corr_true_list = list()
		zscore_true_list = list()
		prob_false_list = list()
		corr_false_list = list()
		zscore_false_list = list()
		
		for row in rows:
			if row["pred_sf"] == row["true_sf"]:
				prob_true_list.append(row["avg_prob"])
				corr_true_list.append(row["corr_coeff"])
				zscore_true_list.append(row["zscore"])
			else:
				#print row["intersection"], self.intersection_limit
				if int(row["intersection"]) > self.intersection_limit:

					prob_false_list.append(row["avg_prob"])
					corr_false_list.append(row["corr_coeff"])
					zscore_false_list.append(row["zscore"])
			
		return prob_false_list, corr_false_list, zscore_false_list

	def plot_SS_functionCorrelation_Analysis(self,list1, list2,c='b', label=""):
		pylab.plot(list1, list2,'o',c=c,label=label,markersize=3)

	def lsqfit(self,list1,list2):
		A = numpy.vstack([list1,numpy.ones(len(list1))]).T
		m,c = numpy.linalg.lstsq(A,list2)[0]
		return m,c

	def plot_by_class(self):
		self.by_class = 'a'
		plist, p_true_list, p_false_list = self.runSSFreqAnalysis()
		#kdrew: use the same bins for all subsequent plots
		a_n,bins = self.plot_SS_Freq_Analysis(plist,c=self.classA_color, marker='o', label="SCOP Class A", bins=self.plot_bins)
		a_n_corr,bins = self.plot_SS_Freq_Analysis(p_true_list,c=self.classA_color, marker='o', label="Class A correct", linetype="--",bins=bins)

		self.by_class = 'b'
		plist, p_true_list, p_false_list = self.runSSFreqAnalysis()
		b_n,bins = self.plot_SS_Freq_Analysis(plist,c=self.classB_color, marker='v', label="SCOP Class B",bins=bins)
		b_n_corr,bins = self.plot_SS_Freq_Analysis(p_true_list,c=self.classB_color, marker='v', label="Class B correct", linetype="--",bins=bins)

		self.by_class = 'c'
		plist, p_true_list, p_false_list = self.runSSFreqAnalysis()
		c_n,bins = self.plot_SS_Freq_Analysis(plist,c=self.classC_color, marker='s', label="SCOP Class C",bins=bins)
		c_n_corr,bins = self.plot_SS_Freq_Analysis(p_true_list,c=self.classC_color, marker='s', label="Class C correct", linetype="--",bins=bins)

		self.by_class = 'd'
		plist, p_true_list, p_false_list = self.runSSFreqAnalysis()
		print len(plist)
		d_n,bins = self.plot_SS_Freq_Analysis(plist,c=self.classD_color, marker='*', label="SCOP Class D",bins=bins)
		d_n_corr,bins = self.plot_SS_Freq_Analysis(p_true_list,c=self.classD_color, marker='*', label="Class D correct", linetype="--",bins=bins)

		#kdrew: o is for other
		self.by_class = 'o'
		plist, p_true_list, p_false_list = self.runSSFreqAnalysis()
		print len(plist)
		o_n,bins = self.plot_SS_Freq_Analysis(plist,c=self.classD_color, marker='*', label="",bins=bins, visible=False)
		o_n_corr,bins = self.plot_SS_Freq_Analysis(p_true_list,c=self.classD_color, marker='*', label="", linetype="--",bins=bins,visible=False)

		pylab.suptitle("Distribution of predicted superfamilies")
		pylab.xlabel("MCM Score")
		pylab.ylabel("Counts")
		leg = pylab.legend(loc=2)
		for t in leg.get_texts():
			t.set_fontsize('small')
		pylab.savefig(self.sincesolved_freq_plot_filename, dpi=150)

		return a_n,a_n_corr,b_n,b_n_corr,c_n,c_n_corr,d_n,d_n_corr,o_n,o_n_corr

	def plot_cor_VS_incor_by_class(self):

		pylab.subplot(221)
		self.by_class = 'a'
		plist, p_true_list, p_false_list = self.runSSFreqAnalysis()
		#kdrew: use the same bins for all subsequent plots
		label_a_correct = r'Class A correct ($\alpha$)'
		a_n_corr,bins = self.plot_SS_Freq_Analysis(p_true_list,c=self.classA_color, marker='o', label=label_a_correct,bins=self.plot_bins)
		a_n_incorr,bins = self.plot_SS_Freq_Analysis(p_false_list,c=self.classA_color, marker='o', label="Class A incorrect", linetype="--", bins=bins, ymax=120)
		leg = pylab.legend(loc=2)
		for t in leg.get_texts():
			t.set_fontsize('small')
		pylab.ylabel("Counts")
		#kdrew: commenting out ratio plotting
		#ax2 = pylab.twinx()
		#pylab.ylabel('ratio')
		#ax2.yaxis.tick_right()
		#ratio = numpy.array(a_n_corr) / numpy.array(a_n)
		#print "ratio: ", ratio

		#l = pylab.plot(bins[1:], ratio, ":", color=self.classA_color,label="Class A ratio",linewidth=1, visible=True)
		#pylab.ylim(ymax=5)
		##a_n_ratio,bins = self.plot_SS_Freq_Analysis(ratio,c=self.classA_color, marker='o', label="Class A ratio", linetype=":",bins=bins, ymax = 4)

		pylab.subplot(223)
		self.by_class = 'b'
		plist, p_true_list, p_false_list = self.runSSFreqAnalysis()
		label_b_correct = r'Class B correct ($\beta$)'
		print "plot correct:"
		b_n_corr,bins = self.plot_SS_Freq_Analysis(p_true_list,c=self.classB_color, marker='v', label=label_b_correct, bins=bins, ymax=30)
		print "plot incorrect:"
		b_n_incorr,bins = self.plot_SS_Freq_Analysis(p_false_list,c=self.classB_color, marker='v', label="Class B incorrect",linetype="--",bins=bins,ymax=30)
		leg = pylab.legend(loc=2)
		for t in leg.get_texts():
			t.set_fontsize('small')
		pylab.ylabel("Counts")
		pylab.xlabel("MCM Score")

		pylab.subplot(224)
		self.by_class = 'c'
		plist, p_true_list, p_false_list = self.runSSFreqAnalysis()
		label_c_correct = r'Class C correct ($\alpha / \beta$)'
		c_n_corr,bins = self.plot_SS_Freq_Analysis(p_true_list,c=self.classC_color, marker='s', label=label_c_correct,bins=bins,ymax=30)
		c_n_incorr,bins = self.plot_SS_Freq_Analysis(p_false_list,c=self.classC_color, marker='s', label="Class C incorrect", linetype="--",bins=bins,ymax=30)
		leg = pylab.legend(loc=2)
		for t in leg.get_texts():
			t.set_fontsize('small')
		pylab.xlabel("MCM Score")

		pylab.subplot(222)
		self.by_class = 'd'
		plist, p_true_list, p_false_list = self.runSSFreqAnalysis()
		print len(plist)
		label_d_correct = r'Class D correct ($\alpha + \beta$)'
		d_n_corr,bins = self.plot_SS_Freq_Analysis(p_true_list,c=self.classD_color, marker='*', label=label_d_correct, bins=bins)
		d_n_incorr,bins = self.plot_SS_Freq_Analysis(p_false_list,c=self.classD_color, marker='*', label="Class D incorrect",linetype="--",bins=bins)
		leg = pylab.legend(loc=2)
		for t in leg.get_texts():
			t.set_fontsize('small')

		#kdrew: o is for other
		self.by_class = 'o'
		plist, p_true_list, p_false_list = self.runSSFreqAnalysis()
		print len(plist)
		o_n_corr,bins = self.plot_SS_Freq_Analysis(p_true_list,c=self.classD_color, marker='*', label="",bins=bins,visible=False)
		o_n_incorr,bins = self.plot_SS_Freq_Analysis(p_false_list,c=self.classD_color, marker='*', label="", linetype="--",bins=bins, visible=False)

		pylab.suptitle("Distribution of predicted superfamilies")
		pylab.savefig(self.sincesolved_freq_plot_filename, dpi=150)

		return a_n_incorr,a_n_corr,b_n_incorr,b_n_corr,c_n_incorr,c_n_corr,d_n_incorr,d_n_corr,o_n_incorr,o_n_corr

if __name__ == "__main__":
	main()


