import MySQLdb
from predictor.structure.structure import *
from domain.domain import *
from metric.log_ratios import *
from metric.mutualInformations import *
from bayes_function_prediction.bayes_function_prediction import *
import predictor.term.func_terms
import domain.mcm_domains

#kdrew: swissprot protein domains
#d = Domain(357441,357441)
#d = Domain(333835,333835)
#d2 = Domain(342401,342401)
#d3 = Domain(333849,333849)

#kdrew: ecoli domains
d = Domain(41086,41086)
print d
d2 = Domain(40974,40974)
print d2
d3 = Domain(37769,37769)
print d3

conn = MySQLdb.connect(host="localhost", user="kdrew", passwd="kdrew_nyu", db="ecoli_benchmark")

#d.load_mcm_structures(conn, scale=0.9)
#print d
#d2.load_mcm_structures(conn, scale=0.9)
#print d2
d3.load_mcm_structures(conn, scale=0.9)
print d3


d3.proc_terms.load_terms(conn,d3.parent_sequence_key, "biological_process")
print d3
d3.loc_terms.load_terms(conn,d3.parent_sequence_key, "cellular_component")
print d3
d3.func_terms.load_terms(conn,d3.parent_sequence_key, "molecular_function")
print d3

