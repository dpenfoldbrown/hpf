

import unittest
import os
import sys
import MySQLdb
import getopt
import ConfigParser
#import shelve
import matplotlib
#kdrew: this is so it creates the figures without accessing the windowing system
matplotlib.use('Agg')
import matplotlib.pyplot
import matplotlib.font_manager as font
import matplotlib.mlab as mlab
import networkx as nx


import pylab
from numpy import arange, array
import numpy

#from hpf.function.metric import *
#from hpf.function.structure import *
#from hpf.function.domain import *
#from hpf.function.bayes import *
#from hpf.function.term import *
#from hpf.function.createdbs.astral import *
#from hpf.go import GODAG,_pydot

from hpf.since_solved.scripts.sincesolved_plots import SinceSolved_Analysis

from matplotlib.pyplot import figure, show
from scipy.stats import gaussian_kde

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Since Solved Error Plots"

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
	sse = SinceSolved_Error(config_file, config_section)

	#kdrew: put all error on same plot but run one at a time
	if sse.onlynewSF_flag and sse.onlyDecoyError_flag and sse.onlyModelError_flag:
		sse.onlynewSF_flag = True
		sse.onlyDecoyError_flag = False
		sse.onlyModelError_flag = False
		histdata_newsf, multiple_histdata = calc_hist_data(ssa,sse)

		sse.onlynewSF_flag = False
		sse.onlyDecoyError_flag = True
		sse.onlyModelError_flag = False
		histdata_decoy, multiple_histdata_decoy = calc_hist_data(ssa,sse)

		sse.onlynewSF_flag = False
		sse.onlyDecoyError_flag = False
		sse.onlyModelError_flag = True
		histdata_model, multiple_histdata_model = calc_hist_data(ssa,sse)

		#kdrew: resetting parameters
		sse.onlynewSF_flag = True
		sse.onlyDecoyError_flag = True
		sse.onlyModelError_flag = True
		

	else:
		histdata, multiple_histdata = calc_hist_data(ssa,sse)

	xlabel = "MCM Score"
	ylabel = ssa.metric1

	fig=figure()
	if ssa.bar_plot_flag:
		pylab.subplot(111)
		xvals = []
		for i in range(len(ssa.plot_bins)-1):
			xvals.append(ssa.plot_bins[i+1]-.045)
		pylab.ylabel("counts")
		pylab.xlabel(xlabel)

		if sse.multiple_histdata_flag:
			pylab.bar(xvals, multiple_histdata, width=.03,color='w')
			
		if sse.onlynewSF_flag and sse.onlyDecoyError_flag and sse.onlyModelError_flag:
			sf_color = "#FFCF79"
			decoy_color = "#E5E4D7"
			model_color = "#2C6700"
			sf_patches = pylab.bar(xvals, histdata_newsf, width=.03, color=sf_color,label="New SF error")
			decoy_patches = pylab.bar(xvals, histdata_decoy, bottom = histdata_newsf, width=.03, color=decoy_color, label="Decoy error")
			#pylab.bar(xvals, histdata_decoy, width=.03, color="r")
			model_patches = pylab.bar(xvals, histdata_model, bottom = map(sum, zip(histdata_newsf,histdata_decoy)), width=.03,color=model_color, label="Classifier error")

			[item.set_hatch("\\") for item in sf_patches]
			[item.set_hatch("/") for item in decoy_patches]
			[item.set_hatch("x") for item in model_patches]

			fp = font.FontProperties(size="x-small")
			#pylab.legend(prop=fp)
			pylab.legend(loc=2)

		else:
			pylab.bar(xvals, histdata, width=.03)

		print "plot_bins:", ssa.plot_bins
		print "xvals:", xvals
		#print "histdata:", histdata
		ax = fig.add_subplot(111)

		if not sse.clean_paper_figure_flag:
			pylab.title(opts,fontsize=10)


		matplotlib.pyplot.savefig(ssa.sincesolved_plot_filename, dpi=150)

def calc_hist_data(ssa,sse):
	data = []
	histdata = []
	multiple_histdata = []

	if ssa.by_mcm_flag:
		for i in range(len(ssa.plot_bins)-1):
			ssa.mcm_value_limit = ssa.plot_bins[i]
			ssa.mcm_value_max = ssa.plot_bins[i+1]
			m1list, m1_true_list, m1_false_list, m2list, m2_true_list, m2_false_list, true_rows, false_rows = ssa.runSSAnalysis()
			if ssa.only_true_flag:
				m1data = m1_true_list	
				m2data = m2_true_list	
			elif ssa.only_false_flag:
				m1data = m1_false_list
				m2data = m2_false_list
			else:
				m1data = m1list
				m2data = m2list

			multiple_histdata.append(len(m1data))

			if sse.onlynewSF_flag:
				m1dataNew = list()
				#for md in m1data:
				for i in range(0,len(m1data)):
					#print md, " new?: ", sse.isNewSF(md)
					#if sse.isNewSF(md):
					if sse.isNewSF(false_rows[i]["true_sf"]):
						print "new sf: ", false_rows[i]["true_sf"]
						m1dataNew.append(false_rows[i]["true_sf"])	

				m1data = m1dataNew

			if sse.onlyContinuous_flag:
				m1dataCont = list()
				#for md in m1data:
				for i in range(0,len(m1data)):
					#kdrew: a little kludgy because I have to find the record in the list of all rows
					#for row in rows:
					#	if row["seqkey"] == md:

					if sse.zscore_cutoff == 0:
						#kdrew: this uses the zscore cutoff to be the decoys zscore to the true sf
						zcutoff = false_rows[i]["zscore"]
					else:
						zcutoff = sse.zscore_cutoff

					if sse.isContinuousSpace(zcutoff, false_rows[i]["true_sf"], false_rows[i]["pred_sf"]):
						print m1data[i], false_rows[i]["true_sf"], false_rows[i]["pred_sf"], "true_newSF:", sse.isNewSF(false_rows[i]["true_sf"])
						m1dataCont.append(m1data[i])

				m1data = m1dataCont

			if sse.onlyModelError_flag:
				m1dataDE = list()
				for i in range(0,len(m1data)):
					#kdrew: m1data should be mcm_zscore and m2data should be zscore to native
					if m1data[i] > m2data[i] and not sse.isNewSF(false_rows[i]["true_sf"]):
						m1dataDE.append(m1data[i])
						print false_rows[i]["seqkey"], m1data[i], m2data[i], false_rows[i]["true_sf"], false_rows[i]["pred_sf"], "true_newSF:", sse.isNewSF(false_rows[i]["true_sf"])

				m1data = m1dataDE

			if sse.onlyDecoyError_flag:
				m1dataDE = list()
				for i in range(0,len(m1data)):
					#kdrew: m1data should be mcm_zscore and m2data should be zscore to native
					if m1data[i] < m2data[i] and not sse.isNewSF(false_rows[i]["true_sf"]):
						m1dataDE.append(m1data[i])
						print false_rows[i]["seqkey"], m1data[i], m2data[i], false_rows[i]["true_sf"], false_rows[i]["pred_sf"], "true_newSF:", sse.isNewSF(false_rows[i]["true_sf"])

				m1data = m1dataDE

			histdata.append(len(m1data))
			if len(m1data) > 0:
				data.append(array(m1data))
			print "Total datapoints: ", sum(histdata)

	return histdata, multiple_histdata


class SinceSolved_Error():
	
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

		self.onlynewSF_flag = config.getboolean(section,'onlynewSF_flag')
		self.onlyContinuous_flag = config.getboolean(section,'onlyContinuous_flag')
		self.onlyDecoyError_flag = config.getboolean(section,'onlyDecoyError_flag')
		self.onlyModelError_flag = config.getboolean(section,'onlyModelError_flag')
		self.multiple_histdata_flag = config.getboolean(section,'multiple_histdata_flag')
		self.zscore_cutoff = config.getfloat(section,'zscore_cutoff')

		self.clean_paper_figure_flag = config.getboolean(section,'clean_paper_figure_flag')


		self.since_solved_conn = MySQLdb.connect(host=self.since_solved_host, user=self.sql_user, passwd=self.sql_password, port=self.sql_port, db=self.since_solved_db)

	def isNewSF(self, sf):
		cursor = self.since_solved_conn.cursor(MySQLdb.cursors.DictCursor)

		query = """
				SELECT * from pdb.new_superfamilies_1_75 where sccs = %s
			"""

		#print query
		cursor.execute(query,(sf, ))
		rows = cursor.fetchall ()
		if len(rows) > 0:
			return True
		else:
			return False

	def isContinuousSpace(self, zscore_cutoff, true_sf, pred_sf):

		zscore_pred2true = self.getZscore(true_sf,pred_sf)

		#print zscore_cutoff, " : ", zscore_pred2true

		#kdrew: this tests to see if the the pred_sf and true_sf are close in structural similarity
		if zscore_cutoff > zscore_pred2true:
			#print "not continuous space"
			return False
		#kdrew: if the pred_sf and the true_sf are closer than the cutoff
		else:
			#print "continuous space"
			return True

		

	def getZscore(self, true_sf, pred_sf):
		cursor = self.since_solved_conn.cursor(MySQLdb.cursors.DictCursor)

		query = """
				SELECT * from pdb.mammoth where experiment_sccs = %s and prediction_sccs = %s
			"""

		#print query
		cursor.execute(query,(true_sf, pred_sf))
		row = cursor.fetchone()
		z_score = -1
		if row == None:
			cursor.execute(query,(pred_sf, true_sf))
			row = cursor.fetchone()
			if row == None:
				return None
			else:
				z_score = row["zscore"]
		else:
			z_score = row["zscore"]

		return z_score


if __name__ == "__main__":
	main()



