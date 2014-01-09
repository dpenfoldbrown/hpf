
import unittest
import os
import MySQLdb

from hpf.function.metric import *
from hpf.function.structure import *
from hpf.function.domain import *
from hpf.function.bayes import *
from hpf.function.term import *

class FunctionPredictionTestCase(unittest.TestCase):
	
	def _metrics(self):
		lr = LogRatios()
		lr.set_metric('MF',None, 1)
		lr.set_metric('MF','P', 10.0)
		lr.set_metric('MF','L', 10.0)
		lr.set_metric('MF','A', 30.0)
		lr.set_metric('MF','a.1.1', 10.0)
		
		mi = MutualInformation()
		#mi.metric_dict[lr._key(('MF',None))] = 0
		mi.set_metric('MF','P', 10.0)
		mi.set_metric('MF','L', 10.0)
		mi.set_metric('MF','Z', 20.0)
		mi.set_metric('MF','A', 30.0)
		mi.set_metric('MF','a.1.1', 10.0)
		mi.set_metric('MF','b.1.1', 20.0)	
		fterms = FuncTerms()
		fterms.append(FuncTerm('MF'))
		return (lr,mi,fterms)
		
	
	def setUp(self):
		self.metrics = self._metrics()
		lr,mi,fterms = self.metrics
		self.bfp = _bfp = BayesFunctionPredictionDB(lr,mi,fterms)

	def testFuncPred(self):
		domain = Domain()
		domain.proc_terms = Terms()
		domain.proc_terms.append(Term('P'))
		domain.loc_terms = Terms()
		domain.loc_terms.append(Term('L'))
		domain.structures = Structures()
		domain.structures.append(Structure('a.1.1',1.0))
		r = self.bfp.predictDomain(domain, -9999999, 9999999)
		assert len(r.predictions) == 1
		assert r.predictions[0].pls_score == 31

		domain = Domain()
		domain.proc_terms = Terms()
		domain.proc_terms.append(Term('P'))
		domain.loc_terms = Terms()
		domain.structures = Structures()
		domain.structures.append(Structure('a.1.1',1.0))		
		r = self.bfp.predictDomain(domain, -9999999, 9999999)
		assert len(r.predictions) == 1
		assert r.predictions[0].pls_score == 10 + 0 + 10 + 1, r.predictions[0].pls_score
		
		domain = Domain()
		domain.proc_terms = Terms()
		domain.proc_terms.append(Term('Z'))
		domain.proc_terms.append(Term('P'))
		domain.loc_terms = Terms()
		domain.structures = Structures()
		domain.structures.append(Structure('a.1.1',1.0))		
		r = self.bfp.predictDomain(domain, -9999999, 9999999)
		assert len(r.predictions) == 1
		assert r.predictions[0].pls_score == LogRatios.DEFAULT_MIN + 0 + 10 + 1, r.predictions[0].pls_score
		
		domain = Domain()
		domain.proc_terms = Terms()
		domain.proc_terms.append(Term('Z'))
		domain.proc_terms.append(Term('P'))
		domain.loc_terms = Terms()
		domain.loc_terms.append(Term('A'))
		domain.loc_terms.append(Term('L'))
		domain.structures = Structures()
		domain.structures.append(Structure('b.1.1',1.0))		
		r = self.bfp.predictDomain(domain, -9999999, 9999999)
		assert len(r.predictions) == 1
		assert r.predictions[0].pls_score == LogRatios.DEFAULT_MIN + 30 + LogRatios.DEFAULT_MIN + 1, r.predictions[0].pls_score
		
		

#class SimpleFunctionPredictionTestCase(unittest.TestCase):
#	def setUp(self):
#		self.lr_metric = LogRatios()
#		self.mi_metric = MutualInformations()
#
#		benchmark_name = "test_benchmark"
#		sql_user="kdrew"
#		sql_password="kdrew_nyu"
#		sql_host="mcpeepants.bio.nyu.edu"
#
#		store_file = "/Users/kdrew/scripts/function_prediction_python/bayes_function_prediction/test_data/test_predictions"
#
#		lr_table_name = "test_log_ratio"
#		mi_table_name = "test_mutual_info"
#
#		self.hpf_conn = MySQLdb.connect(host=sql_host, user=sql_user, passwd=sql_password, db="hpf")
#		self.mygo_conn = MySQLdb.connect(host=sql_host, user=sql_user, passwd=sql_password, db="mygo")
#		self.data_conn = MySQLdb.connect(host=sql_host, user=sql_user, passwd=sql_password, db=benchmark_name)
#
#		self.lr_metric.load_metric(self.hpf_conn, lr_table_name)
#
#		self.mi_metric.load_metric(self.hpf_conn,mi_table_name)
#
#		fterms = predictor.term.func_terms.FuncTerms(self.mygo_conn)
#
#		self.domains = domain.mcm_domains.MCMDomains(self.data_conn, scale=0.8)
#
#		bfp = BayesFunctionPredictionFS(self.lr_metric,self.mi_metric, fterms)
#		bfp.predictBenchmarkDomains(self.domains, eagerStore=0)
#		bfp.storeBenchmarkPredictions(self.domains,filename=store_file)
#
#
#	def testFuncPred(self):
#		domain = self.domains[0]
#		print domain.domain_sequence_key, ", ", len(domain.predictions)
#		prediction = domain.predictions.getByAcc("GO:0015002")
#
#		assert 7457 == len(domain.predictions), "wrong number of predictions"
#		assert round(13.43857,3) == round(prediction.s_score,3), "wrong structure score"
#		assert round(1.0986,3) == round(prediction.p_score,3), "wrong process score"

if __name__ == "__main__":
            unittest.main()
    


