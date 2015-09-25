
import unittest
import os
import MySQLdb
import getopt
import ConfigParser
import shelve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import networkx as nx
#kdrew: this is so it creates the figures without accessing the windowing system


import pylab
from numpy import arange

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
	print "Welcome to SCOP Superfamily Function Analysis"

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

	sfa = SF_Func_Analysis(config_file, config_section)
	if sfa.do_sf_prob_analysis:
		ls, ps = sfa.runSFProbAnalysis()
		sfa.plotSFProbAnalysis(ls, ps)
	if sfa.do_sf_func_analysis:
		sfa.runSFFuncAnalysis()
	if sfa.do_sf_func_analysisOne:
		print "Running One"
		sfa.runSFFuncOne(sf=sfa.superfamily)
	

class SF_Func_Analysis():
	
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
		self.astral_db = config.get(section, 'astral_db')
		self.astral_host = config.get(section, 'astral_host')
		self.astral_table = config.get(section, 'astral_table')

		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')
		self.prob_table_name = config.get(section, 'probability_table_name')

		self.noLoad = config.getboolean(section,'noLoad')
		self.probability_shelve = config.get(section,'probability_shelve')

		self.do_sf_prob_analysis = config.getboolean(section,'do_sf_prob_analysis')
		self.sf_prob_plot_filename = config.get(section,'sf_prob_plot_filename')

		self.do_sf_func_analysis = config.getboolean(section,'do_sf_func_analysis')
		self.do_sf_func_analysisOne = config.getboolean(section,'do_sf_func_analysisOne')
		self.sf_func_plot_filename = config.get(section,'sf_func_plot_filename')
		self.superfamily = config.get(section,'superfamily')
		

		dict = shelve.open(self.probability_shelve)
		
		self.prob_metric = Probability(dict=dict)

		self.mygo_conn = MySQLdb.connect(host=self.mygo_host, user=self.sql_user, passwd=self.sql_password, db=self.mygo_db)

		if not self.noLoad:
			self.metric_conn = MySQLdb.connect(host=self.metric_host, user=self.sql_user, passwd=self.sql_password, db=self.metric_db)
			self.prob_metric.load_metric(self.metric_conn, self.prob_table_name)

		self.astl = Astral(self.sql_user,self.sql_password, self.astral_host, db_name=self.astral_db, table=self.astral_table)
		self.astl.connect()

	#kdrew: creates piecharts
	def runSFFuncAnalysis(self):
		#kdrew: get global probabilities for all superfamilies
		for sf in self.astl.get_all_superfamilies()[1089:]: #temporarily skipping ahead
			self.runSFFuncOne(sf)

	def runSFFuncOne(self, sf):
		#print sf, ": ", self.prob_metric.get_metric("all",sf)
		sf_prob_met = self.prob_metric.get_all_metric(sf)
		go_labels = list()
		go_probs = list()
		go_acc = list()
		all_labels = list()
		for key in sf_prob_met:
			#kdrew: only look at molecular_function terms and under a certain background prob
			if is_type(key[0],"molecular_function",self.mygo_conn):
				all_labels.append(key[0]+"\n"+get_name(key[0], self.mygo_conn))
				if 0.1 >= self.prob_metric.get_metric("all",key[0]):
					#print sf,key[0],": ",sf_prob_met.get_metric(key[0])
					go_labels.append(key[0]+"\n"+get_name(key[0], self.mygo_conn))
					go_probs.append(sf_prob_met.get_metric(key[0]))
					go_acc.append(key[0])
				
		self.plotSFFuncAnalysis(sf, go_labels,go_probs)
		self.plotSFGOTree(sf,[Term(label.split("\n")[0],label.split("\n")[1]) for label in all_labels])

	def plotSFGOTree(self, sf, terms):
		dag = GODAG(nodes=terms)
		#print "Loading edges"
		dag.load_edges(self.mygo_conn)
		node_color = {}
		for node in dag.nodes():
			node_color[_pydot(node)] = self.prob_metric.get_metric("all",node) 
		node_size = {}
		for node in dag.nodes():
			node_size[_pydot(node)] = 50 + self.prob_metric.get_metric(sf, node) * 1000 
		#print dag.nodes()
		graph = dag.graph()
		#print graph.nodes()
		labels = {}
		for node in graph.nodes():
			labels[node] = "GO:"+node.acc+"\n"+node.name
			
		# Now add the ruler
		rang = range(11)
		rang.reverse()
		base_prob =  [n*0.1 for n in rang]
		base_lab = ["%i" % int(n*10) for n in base_prob]
		for i,l in enumerate(base_lab):
			node_color[l] = base_prob[i]
			node_size[l] = 50 + base_prob[i]*1000
			labels[l] = "0."+l if l != "10" else "1"
			graph.add_node(l)
			if i < 10:
				graph.add_edge(l,base_lab[i+1])
		#print graph.nodes()
		
		pos = nx.pydot_layout(graph, prog="dot")
		matplotlib.pyplot.figure(figsize=(10, 10))
		nx.draw(graph,
			    font_size=8,
			    labels=labels,
			    pos=pos,
			    cmap="autumn",
			    node_color=[node_color[n] for n in graph.nodes()],
			    node_size=[node_size[n] for n in graph.nodes()],
			    with_labels=True)
		f_parts = self.sf_func_plot_filename.rpartition('.')
		filename = f_parts[0] + "_" + sf + "_go_" + "." + f_parts[2]
		matplotlib.pyplot.savefig(filename)

	def plotSFFuncAnalysis(self,sf,go_labels,go_probs):
		pylab.figure(figsize=(8,8))
		ax = pylab.axes([0.3, 0.5, 0.4, 0.4])
		patches, texts = pylab.pie(go_probs,labels=go_labels)
		labels=go_labels
		pylab.title("Function Frequency of SCOP superfamilies: "+ sf)
		#kdrew: a little filename manipulation
		f_parts = self.sf_func_plot_filename.rpartition('.')
		pylab.savefig(f_parts[0]+"_"+sf+"."+f_parts[2])

	#kdrew: creates a barchart of frequencies of superfamilies across all sequences used to train probability tables
	def runSFProbAnalysis(self):
		sf_labels = list()
		sf_probs = list()
		#kdrew: get global probabilities for all superfamilies
		for sf in self.astl.get_all_superfamilies():
			#print sf, ": ", self.prob_metric.get_metric("all",sf)
			sf_labels.append(sf)
			sf_probs.append(self.prob_metric.get_metric("all",sf))

		return (sf_labels,sf_probs)

	def plotSFProbAnalysis(self, sf_labels, sf_probs):
		rects = pylab.bar(arange(len(sf_labels)),sf_probs)
		pylab.ylim(ymax=.1, ymin=0.0)
		pylab.title("Frequency of SCOP superfamilies")
		pylab.xlabel("SCOP superfamilies")
		pylab.ylabel("Frequency of sequences (%)")
		self.autolabel(rects,sf_labels)
		pylab.savefig(self.sf_prob_plot_filename)
			
	def autolabel(self, rects, labels):
		i = 0
		for rect in rects:
			i = i+1
			height = rect.get_height()
			if height > .01:
				pylab.text(rect.get_x()+rect.get_width()/2.,1.05*height, labels[i-1],ha='center',va='bottom',fontsize=6, rotation=45)


if __name__ == "__main__":
	main()
    


