
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


import pylab
from numpy import arange, array
import numpy

from hpf.function.metric import *
from hpf.function.structure import *
from hpf.function.domain import *
from hpf.function.bayes import *
from hpf.function.term import *
from hpf.function.createdbs.astral import *
from hpf.go import GODAG,_pydot

from matplotlib.pyplot import figure, show
from scipy.stats import gaussian_kde

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Since Solved Plots"

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

	data = []
	histdata = []
	if ssa.by_class_flag:
		ssa.by_class = 'a'
		plist, p_true_list, p_false_list, mlist, m_true_list, m_false_list, true_rows, false_rows = ssa.runSSAnalysis()
		if ssa.violin_plot_flag:
			data.append(array(mlist))

		if ssa.scatter_plot_flag:
			if not ssa.false_only_flag:
				pylab.scatter(p_true_list, m_true_list, c=ssa.classA_color, marker='o')
				pylab.scatter(p_false_list, m_false_list, c=ssa.classA_color, marker='x')
			else:
				pylab.scatter(p_false_list, m_false_list, c=ssa.classA_color, marker='o')

		ssa.by_class = 'b'
		plist, p_true_list, p_false_list, mlist, m_true_list, m_false_list,rows = ssa.runSSAnalysis()
		if ssa.violin_plot_flag:
			data.append(array(mlist))
		if ssa.scatter_plot_flag:
			if not ssa.false_only_flag:
				pylab.scatter(p_true_list, m_true_list, c=ssa.classB_color, marker='o')
				pylab.scatter(p_false_list, m_false_list, c=ssa.classB_color, marker='x')
			else:
				pylab.scatter(p_false_list, m_false_list, c=ssa.classB_color, marker='o')

		ssa.by_class = 'c'
		plist, p_true_list, p_false_list, mlist, m_true_list, m_false_list,true_rows, false_rows = ssa.runSSAnalysis()
		if ssa.violin_plot_flag:
			data.append(array(mlist))
		if ssa.scatter_plot_flag:
			if not ssa.false_only_flag:
				pylab.scatter(p_true_list, m_true_list, c=ssa.classC_color, marker='o')
				pylab.scatter(p_false_list, m_false_list, c=ssa.classC_color, marker='x')
			else:
				pylab.scatter(p_false_list, m_false_list, c=ssa.classC_color, marker='o')

		ssa.by_class = 'd'
		plist, p_true_list, p_false_list, mlist, m_true_list, m_false_list,true_rows, false_rows = ssa.runSSAnalysis()
		print len(plist)
		if ssa.violin_plot_flag:
			data.append(array(mlist))
		if ssa.scatter_plot_flag:
			if not ssa.false_only_flag:
				pylab.scatter(p_true_list, m_true_list, c=ssa.classD_color, marker='o')
				pylab.scatter(p_false_list, m_false_list, c=ssa.classD_color, marker='x')
			else:
				pylab.scatter(p_false_list, m_false_list, c=ssa.classD_color, marker='o')
	elif ssa.by_correct_flag:
		m1list, m1_true_list, m1_false_list, m2list, m2_true_list, m2_false_list,true_rows, false_rows = ssa.runSSAnalysis()
		if ssa.violin_plot_flag:
			data.append(array(m1_true_list))
			data.append(array(m1_false_list))

		if ssa.scatter_plot_flag:
				pylab.scatter(m1_true_list, m2_true_list, c='b')
				pylab.scatter(m1_false_list, m2_false_list, c='r')

	elif ssa.by_mcm_flag:
		running_sum = 0
		running_num = 0
		for i in range(len(ssa.plot_bins)-1):
			ssa.mcm_value_limit = ssa.plot_bins[i]
			ssa.mcm_value_max = ssa.plot_bins[i+1]
			m1list, m1_true_list, m1_false_list, m2list, m2_true_list, m2_false_list,true_rows, false_rows = ssa.runSSAnalysis()
			if ssa.only_true_flag:
				m1data = m1_true_list	
			elif ssa.only_false_flag:
				m1data = m1_false_list
			else:
				m1data = m1list
			running_num+= len(m1data)
			running_sum+= sum(m1data)
			if ssa.violin_plot_flag:
				histdata.append(len(m1data))
				if len(m1data) > 0:
					data.append(array(m1data))
				print "Total datapoints: ", sum(histdata)
				xlabel = "MCM Score"
				ylabel = ssa.metric1


			if ssa.scatter_plot_flag:
				pylab.scatter(m1_true_list, m2_true_list, c='b')
				pylab.scatter(m1_false_list, m2_false_list, c='r')

		print "average: ", running_sum*1.0/running_num

	else:
		plist, p_true_list, p_false_list, mlist, m_true_list, m_false_list,true_rows, false_rows = ssa.runSSAnalysis()
		pylab.scatter(p_true_list, m_true_list, c='b')
		pylab.scatter(p_false_list, m_false_list, c='r')

	if ssa.violin_plot_flag:
		fig=figure()
		if ssa.bar_plot_flag:
			pylab.subplot(413)
			xvals = []
			for i in range(len(ssa.plot_bins)-1):
				#xvals.append(ssa.plot_bins[i+1]-.045)
				xvals.append(ssa.plot_bins[i+1]-.065)
			pylab.ylabel("counts")
			pylab.xlabel(xlabel)
			pylab.bar(xvals, histdata, width=.03)
			print "plot_bins:", ssa.plot_bins
			print "xvals:", xvals
			print "histdata:", histdata
			ax = fig.add_subplot(211)
		else:
			ax = fig.add_subplot(111)

		if not ssa.paper:
			pylab.title(opts,fontsize=10)
		pylab.ylabel(ylabel)
		ssa.violin_plot(ax,data,pos=range(len(data)),bp=True)
		blanks = list("" for i in range(len(data)))
		pylab.xticks(range(len(data)), blanks)
		matplotlib.pyplot.savefig(ssa.sincesolved_plot_filename, dpi=150)
	else:
		pylab.title(opts,fontsize=10)
		pylab.suptitle("Distribution of predicted superfamilies")
		pylab.xlabel(ssa.metric1)
		pylab.ylabel(ssa.metric2)
		#leg = pylab.legend(loc=2)
		#for t in leg.get_texts():
		#	t.set_fontsize('small')
		pylab.savefig(ssa.sincesolved_plot_filename, dpi=300)


class SinceSolved_Analysis():
	
	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		self.since_solved_db = config.get(section, 'since_solved_db')
		self.since_solved_host = config.get(section, 'since_solved_host')
		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')
		self.sql_port= config.getint(section, 'sql_port')

		self.sincesolved_plot_filename = config.get(section,'sincesolved_plot_filename')

		self.mcm_value_limit = config.getfloat(section,'mcm_value_limit')
		self.mcm_value_max = config.getfloat(section,'mcm_value_max')
		self.scop_level = config.getint(section,'scop_level')
		self.plot_bins_str = config.get(section,'plot_bins')
		if self.plot_bins_str == "":
			self.plot_bins = 10
		else:
			self.plot_bins = map(float,list(self.plot_bins_str.split(',')))
		self.by_class_flag = config.getboolean(section,'by_class_flag')
		self.by_correct_flag = config.getboolean(section,'by_correct_flag')
		self.by_mcm_flag = config.getboolean(section,'by_mcm_flag')
		self.only_true_flag = config.getboolean(section,'only_true_flag')
		self.only_false_flag = config.getboolean(section,'only_false_flag')
		self.group_sf_flag = config.getboolean(section,'group_sf_flag')
		self.by_class = ""
		self.date_cutoff = config.get(section,'date_cutoff')
		self.exclude_clusters_str = config.get(section,'exclude_clusters')
		self.exclude_clusters = []
		if self.exclude_clusters_str != "":
			self.exclude_clusters = map(int, list(self.exclude_clusters_str.split(',')))

		self.metric1 = config.get(section,'metric1')
		self.metric2 = config.get(section,'metric2')
		self.false_only_flag = config.getboolean(section,'false_only_flag')

		self.violin_plot_flag = config.getboolean(section,'violin_plot_flag')
		self.bar_plot_flag = config.getboolean(section,'bar_plot_flag')
		self.scatter_plot_flag = config.getboolean(section,'scatter_plot_flag')

		self.sample_denovo_flag =  config.getboolean(section,'sample_denovo_flag')
		self.sample_denovo = config.getint(section,'sample_denovo')
		self.sample_denovo_rand_seed = config.getint(section,'sample_denovo_rand_seed')

		self.redux = config.getboolean(section,'redux')
		self.paper = config.getboolean(section,'paper')

		self.classA_color = "blue"
		self.classB_color = "orange"
		self.classC_color = "red"
		self.classD_color = "green"

		self.since_solved_conn = MySQLdb.connect(host=self.since_solved_host, user=self.sql_user, passwd=self.sql_password, port=self.sql_port, db=self.since_solved_db)

	def runSSAnalysis(self):
		cursor = self.since_solved_conn.cursor(MySQLdb.cursors.DictCursor)
		if self.by_class != "":
			by_class_query = " and SUBSTRING_INDEX(mc.experiment_sccs,'.',1) = '"+self.by_class+"'"
		else:
			by_class_query = ""

		if self.group_sf_flag:
			group_sf_query = " GROUP BY SUBSTRING_INDEX(mc.experiment_sccs,'.',%s) " % (self.scop_level,)
		else:
			group_sf_query = ""

		query = """
				SELECT mc.cluster_id AS seqkey, SUBSTRING_INDEX(mc.experiment_sccs,'.',%s) AS pred_sf, SUBSTRING_INDEX(ss.sccs,'.',%s) AS true_sf, mc.probability as top_prob,
					m.*, mc.zscore as mcm_zscore, mc.convergence, mc.prediction_contact_order, length(s.sequence) as seqlen, length(s.sequence)*(m.end_psi/100) as resAli, max_table.avg_mcm_zscore, max_table.min_mcm_zscore
				FROM since_solved.mcm_cluster mc, since_solved.since_solved ss, 
					(SELECT cluster_id, MAX(probability) as max_prob, AVG(zscore) as avg_mcm_zscore, MIN(zscore) as min_mcm_zscore from since_solved.mcm_cluster group by cluster_id) as max_table ,
					since_solved.mammoth as m, hpf.sequence as s
				WHERE mc.cluster_id = ss.cluster_id and max_table.cluster_id = mc.cluster_id and mc.probability = max_table.max_prob 
					AND m.prediction = mc.structure_key
					AND mc.sequence_key = s.id
					AND ss.ascession_date >= DATE(%s)
					AND mc.probability > %s AND mc.probability <= %s
			""" + by_class_query + group_sf_query + """ 
				ORDER BY top_prob
			"""
		if self.redux:
			if self.paper:
				query = """ 
						SELECT mc.*, cluster_id as seqkey, 
							m.*, length(s.sequence) as seqlen, length(s.sequence)*(m.end_psi/100) as resAli
						FROM since_solved.p_mcm_redux as mc, since_solved.mammoth_redux as m, hpf.sequence as s
						WHERE m.prediction = mc.structure_key AND mc.sequence_key = s.id 
							AND mc.top_prob > %s AND mc.top_prob <= %s
						""" + by_class_query +""" 
						ORDER BY cluster_id 
					"""

		if self.sample_denovo_flag:
			query = """
				SELECT * FROM 
					(SELECT md.sequence_key AS seqkey, SUBSTRING_INDEX(md.experiment_sccs,'.',%s) AS pred_sf, md.probability AS top_prob, 
						md.zscore as mcm_zscore, md.convergence, md.prediction_contact_order 
					FROM mcmData AS md, (SELECT sequence_key, MAX(probability) AS max_prob FROM mcmData GROUP BY sequence_key) AS max_table, 
						hpf.domain AS d, hpf.protein AS p, p_experiments AS e 
					WHERE md.sequence_key = d.domain_sequence_key AND d.parent_sequence_key = p.sequence_key AND p.experiment_key = e.id 
						AND max_table.sequence_key = md.sequence_key AND max_table.max_prob = md.probability 
					ORDER BY RAND(%s) LIMIT %s) as tmp
					WHERE top_prob > %s AND top_prob <= %s
				"""
			print query
			cursor.execute(query,(self.scop_level,self.sample_denovo_rand_seed, self.sample_denovo, self.mcm_value_limit,self.mcm_value_max ))

		else:
			print query
			if self.paper:
				cursor.execute(query,(self.mcm_value_limit,self.mcm_value_max))
			else:
				cursor.execute(query,(self.scop_level,self.scop_level,self.date_cutoff,self.mcm_value_limit,self.mcm_value_max))

		rows = cursor.fetchall ()
	
		metric1_true_list = list()
		metric1_false_list = list()
		metric1_list = list()
		metric2_true_list = list()
		metric2_false_list = list()
		metric2_list = list()
		true_rows = list()
		false_rows = list()
		
		for row in rows:
			if self.sample_denovo_flag:
				#print row["seqkey"], row["pred_sf"], row["top_prob"], row["mcm_zscore"], row["convergence"], row["prediction_contact_order"]
				metric1_list.append(row[self.metric1])
				metric2_list.append(row[self.metric2])
			else:
				if int(row["seqkey"]) not in self.exclude_clusters:
				#	print row["seqkey"], row["pred_sf"], row["true_sf"], row["top_prob"], row["zscore"], row["ini_rms"], row["convergence"], row["prediction_contact_order"]

					if row["pred_sf"] == row["true_sf"]:
						metric1_true_list.append(row[self.metric1])
						metric2_true_list.append(row[self.metric2])
						true_rows.append(row)
					else:
						metric1_false_list.append(row[self.metric1])
						metric2_false_list.append(row[self.metric2])
						false_rows.append(row)

					metric1_list.append(row[self.metric1])
					metric2_list.append(row[self.metric2])

		print len(metric1_list)
		return metric1_list, metric1_true_list, metric1_false_list, metric2_list, metric2_true_list, metric2_false_list, true_rows, false_rows

	def violin_plot(self,ax,data,pos, bp=False):
	    '''
	    create violin plots on an axis
	    '''
	    dist = max(pos)-min(pos)
	    w = min(0.15*max(dist,1.0),0.5)
	    for d,p in zip(data,pos):
		print d,p
		k = gaussian_kde(d) #calculates the kernel density
		m = k.dataset.min() #lower bound of violin
		M = k.dataset.max() #upper bound of violin
		x = arange(m,M,(M-m)/100.) # support for violin
		v = k.evaluate(x) #violin profile (density curve)
		v = v/v.max()*w #scaling the violin to the available space
		ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
		ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
	    if bp:
		ax.boxplot(data,notch=1,positions=pos,vert=1)


if __name__ == "__main__":
	main()


