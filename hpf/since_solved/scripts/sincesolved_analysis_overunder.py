
import unittest
import os
import sys
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

from hpf.since_solved.scripts.sincesolved_analysis import SinceSolved_Analysis

import pylab
from numpy import arange
import numpy
from matplotlib.pyplot import figure, show

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Since Solved Analysis Over Under"

	config_file = "config.cfg"
	config_section = "Default"

	try:
		print "before getopt"
		opts, args = getopt.getopt(sys.argv[1:], "hc:s:S:", ["help", "config=", "section=", "section2="]) 
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
		elif opt in ("-s","--section2"):
			config_section2 = arg
		else:
			assert False, "unhandled option"

	ssaUnder = SinceSolved_Analysis(config_file, config_section)
	if ssaUnder.do_sincesolved_mammoth_analysis:
		ssaUnder.runSSMammothAnalysis()
	if ssaUnder.do_sincesolved_freq_analysis:
		if ssaUnder.by_class_flag:
			#a_n,a_n_corr,b_n,b_n_corr,c_n,c_n_corr,d_n,d_n_corr = ssaUnder.plot_by_class()
			under_results = ssaUnder.plot_by_class()

	ssaOver = SinceSolved_Analysis(config_file, config_section2)
	if ssaOver.do_sincesolved_mammoth_analysis:
		ssaOver.runSSMammothAnalysis()
	if ssaOver.do_sincesolved_freq_analysis:
		if ssaOver.by_class_flag:
			#a_n,a_n_corr,b_n,b_n_corr,c_n,c_n_corr,d_n,d_n_corr = ssaOver.plot_by_class()
			over_results = ssaOver.plot_by_class()

	print under_results[1][-1:], under_results[0][-1:]
	print under_results[1][-1:]*1.0/ under_results[0][-1:]
	print over_results[1][-1:], over_results[0][-1:]
	print over_results[1][-1:]*1.0/ over_results[0][-1:]

	data = []
	data2 = []
	for i in range(0,len(under_results),2):
		data.append(under_results[i+1][-1:][0]*1.0/ under_results[i][-1:][0])
		data2.append(over_results[i+1][-1:][0]*1.0/ over_results[i][-1:][0])
		
	print data
	print data2

	fig=figure()
	pylab.subplot(111)
	xvals = range(2,9,2)
	print xvals
	width = 0.35
	pylab.bar(xvals,data2,width=width)
	pylab.xticks(xvals, ('A', 'B', 'C', 'D') )
	xvals = range(3,10,2)
	print xvals
	pylab.bar(xvals,data,width=width,color="g")
	pylab.ylabel("Percent correct")
	ax = fig.add_subplot(111)

	pylab.title(opts,fontsize=10)

	matplotlib.pyplot.savefig(ssaUnder.sincesolved_freq_plot_filename+"over_under.png", dpi=300)

if __name__ == "__main__":
	main()


