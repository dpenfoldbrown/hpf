import MySQLdb
from predictor.structure.structure import *
from domain.domain import *
from metric.log_ratios import *
from metric.mutualInformations import *
from bayes_function_prediction import *
import predictor.term.func_terms
import domain.mcm_domains
from domain.mcm_domain import *

#kdrew: swissprot protein domains
#d = Domain(357441,357441)
#d = Domain(333835,333835)
#d2 = Domain(342401,342401)
#d3 = Domain(333849,333849)

#kdrew: ecoli domains
d = MCMDomain(41086,41086)
print d
d2 = MCMDomain(40974,40974)
print d2
d3 = MCMDomain(37769,37769)
print d3

#conn = MySQLdb.connect(host="localhost", user="kdrew", passwd="kdrew_nyu", db="swissprot_benchmark")
conn = MySQLdb.connect(host="localhost", user="kdrew", passwd="kdrew_nyu", db="ecoli_benchmark")

#d.load_structures(conn, scale=0.9)
#print d
#d.structures.print_structures()
#d2.load_structures(conn, scale=0.9)
#print d2
#d3.load_structures(conn, scale=0.9)
#print d3

#print "compare:",compare(d.structures.structs[0],d.structures.structs[1])


#print d3.proc_terms.evidence_codes
d3.proc_terms.load_terms(conn,d3.parent_sequence_key, "biological_process")
d3.proc_terms.print_terms()
d3.loc_terms.load_terms(conn,d3.parent_sequence_key, "cellular_component")
d3.loc_terms.print_terms()
d3.func_terms.load_terms(conn,d3.parent_sequence_key, "molecular_function")
d3.func_terms.print_terms()


#hpf_conn = MySQLdb.connect(host="localhost", user="kdrew", passwd="kdrew_nyu", db="hpf")
#lr = LogRatios()
#lr.load_function_logRatio(hpf_conn)
#lr.load_superfamily_logRatio(hpf_conn,"sfGmf_LR_notEcoli_IEA8")
#lr.load_process_logRatio(hpf_conn,"bpGmf_LR_notEcoli_IEA8")
#lr.load_localization_logRatio(hpf_conn,"ccGmf_LR_notEcoli_IEA8")

#mi = MutualInformations()
#mi.load_process_mutual_information(hpf_conn,"bpMImf_notEcoli_IEA8")
#mi.load_localization_mutual_information(hpf_conn,"ccMImf_notEcoli_IEA8")


#for func in d3.func_terms:
	#print "Featured term"
#	featured_proc = d3.proc_terms.getFeaturedTerm(func,mi)
#	mf_bp_lr = lr.get_logratio(func,featured_proc)
#	mf_lr = lr.get_logratio(func)
	#func.print_term()
	#print featured_proc
	#print "mf_bp_lr: ", mf_bp_lr, " mf_lr: ",mf_lr, " mf_bp_lr+mf_lr", mf_bp_lr+mf_lr
	#for proc in d3.proc_terms:
		#func.print_term()
		#proc.print_term()


#fterms = predictor.term.func_terms.FuncTerms()
#bfp = BayesFunctionPrediction(lr,mi, fterms)
#fterms.load_terms(hpf_conn)

#bfp.computeDomain(d)
#bfp.computeDomain(d2)
#bfp.computeDomain(d3)

#ecoli_domains = domain.mcm_domains.MCMDomains("ecoli_benchmark", conn)

#results_conn = MySQLdb.connect(host="localhost", user="kdrew", passwd="kdrew_nyu", db="hpf_results")
#cursor = results_conn.cursor(MySQLdb.cursors.DictCursor)
#bfp.createResultsTable(cursor,"temp_bayes_results")
#bfp.uploadPredictions(cursor, d2,"temp_bayes_results")


