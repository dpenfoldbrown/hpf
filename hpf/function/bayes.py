import MySQLdb

import sys
from metric import *
from prediction import *
from collections import defaultdict

class BayesFunctionPrediction:

	__SCORE_TYPE = "llr"

	def __init__(self, log_ratio, mutual_info, global_functions):
		self.log_ratio = log_ratio
		self.mutual_info = mutual_info
		self.global_functions = global_functions

	def predictBenchmarkDomains(self, domains):
		for domain in domains:
			self.predictDomain(domain)

	#kdrew: mutates domain by adding new predictions
	def predictDomain(self, domain, lower_bound=-12, upper_bound=0):	
		for func in self.global_functions:
			pred = Prediction(func,self.__SCORE_TYPE)
			s_score_raw = self.computeStruct(domain.structures,func)
			p_score_raw, proc_feature = self.computeTerms(func,domain.proc_terms)
			l_score_raw, loc_feature = self.computeTerms(func,domain.loc_terms)
			f_score_raw, no_feature = self.computeTerms(func)
			assert(f_score_raw != 0)

			pred.update_s_score(f_score_raw + s_score_raw)
			pred.update_p_score(f_score_raw + p_score_raw)
			pred.update_l_score(f_score_raw + l_score_raw)
			pred.update_ps_score(f_score_raw + s_score_raw + p_score_raw)
			pred.update_ls_score(f_score_raw + s_score_raw + l_score_raw)
			pred.update_pl_score(f_score_raw + l_score_raw + p_score_raw)
			pred.update_pls_score(f_score_raw + s_score_raw + p_score_raw + l_score_raw)
			pred.update_base_score(f_score_raw)

			pred.update_proc_feature(proc_feature)
			pred.update_loc_feature(loc_feature)

			#if s_score_raw > -16:
			#	print s_score_raw,p_score_raw,l_score_raw,f_score_raw,pred.pls_score

			
			#if not all([score == LogRatios.DEFAULT_MIN for score in (f_score_raw,p_score_raw,l_score_raw)]):
			#	print pred.pls_score, (f_score_raw,p_score_raw,l_score_raw)
			#if s_score_raw > LogRatios.DEFAULT_MIN + 1:
			#	print domain.structures, s_score_raw, s_score_raw-LogRatios.DEFAULT_MIN
					
			
			#print pred.pls_score
			if pred.pls_score >= lower_bound and pred.base_score < upper_bound:
				domain.append_prediction(pred)
		#print "Predicted",len(domain.predictions)
		return domain

	def computeStruct(self,structures,func):
		s_running_sum = 0
		for struct in structures:
			# Sum the likelihood each structure gives to the function
			# This assumes the structure probabilities sum to 1, are pre-scaled
			s_running_sum += struct.get_prob()*self.log_ratio.get_metric(func,struct)
			#print "struct.get_prob():",struct.get_prob(),"struct.superfamily:",struct.superfamily,"metric:",self.metric.get_metric(struct,func)
		return s_running_sum

	def computeTerms(self,func,terms=None):
		if None == terms:
			# This is to get the base llr for the function
			return self.log_ratio.get_metric(func,None), None
		if len(terms)==0:
			# This is if there aren't any known GO terms.
			return 0,None
		featured_term = terms.getFeaturedTerm(func,self.mutual_info)
		if None != featured_term:
			return self.log_ratio.get_metric(func,featured_term), featured_term
		raise Exception("Shouldn't ever reach here")
				
	def formatPredictions(self, domain):
		"""Returns the best pls prediction for each MF as a list of tuples"""
		best_mf = defaultdict(lambda: -sys.maxint+1)
		mf = {}
		# Search for the best prediction per molecular function
		for pred in domain.predictions:
			func_acc = pred.function.get_id()
			if pred.pls_score > best_mf[func_acc]:
				mf[func_acc] = pred
				best_mf[func_acc] = pred.pls_score
		mf_acc = set(mf.keys())
		#print len(mf_acc)," unique functions predicted"
		
		insert_list = []
		#for pred in mf.values()
		for pred in mf.values():
			mf_acc.add(pred.function.get_id())
			p_acc = pred.proc_feature.get_id() if pred.proc_feature != None else None
			l_acc = pred.loc_feature.get_id() if pred.loc_feature != None else None
			sccs = ",".join([structure.get_id() for structure in domain.structures])
			insert_list.append((domain.parent_sequence_key, domain.domain_sequence_key, pred.function.get_id(), pred.function.name, pred.p_score, pred.l_score, pred.s_score, pred.pl_score, pred.ps_score, pred.ls_score, pred.pls_score, pred.base_score, domain.type, p_acc, l_acc, sccs))
		return insert_list

class BayesFunctionPredictionFS(BayesFunctionPrediction):

	__STORAGE_FILENAME = "bayes_llr_results"

	def predictBenchmarkDomains(self, domains, eagerStore=1, filename=__STORAGE_FILENAME):
		if eagerStore:
			fileHandle = open(filename,"a")
		for domain in domains:
			self.predictDomain(domain)
			if eagerStore:
				self.storePredictions(fileHandle, domain)
				domain.predictions=Predictions(uploaded=1)
		if eagerStore:
			fileHandle.close()

	def storeBenchmarkPredictions(self, domains, filename):
		fileHandle = open(filename,"a")
		for domain in domains:
			self.storePredictions(fileHandle, domain)
		fileHandle.close()
		
	def storePredictions(self,fileHandle, domain):
		insert_string = self.formatPredictions(domain)
		fileHandle.write(insert_string)

	def formatPredictions(self, domain):
		insert_string = ""
		for pred in domain.predictions:
			insert_string = insert_string+ str(domain.parent_sequence_key)+"\t"+str(domain.domain_sequence_key)+"\t"+str(pred.function.get_id())+"\t"+ str(pred.function.name)+"\t"+ str(pred.p_score)+"\t"+ str(pred.l_score )+"\t"+ str(pred.s_score)+"\t"+ str(pred.pl_score)+"\t"+ str(pred.ps_score)+"\t"+ str(pred.ls_score)+"\t"+ str(pred.pls_score)+"\t"+ str(pred.base_score)+"\n"

		return insert_string
	
class BayesFunctionPredictionDB(BayesFunctionPrediction):

	__TEMPLATE_TABLE_NAME = "bayes_llr_results_template"

	def predictBenchmarkDomains(self, domains, eagerUpload_conn=None, eager_table_name=None):
		if None != eagerUpload_conn:
			cursor = eagerUpload_conn.cursor(MySQLdb.cursors.DictCursor)
			self.createResultsTable(cursor,eager_table_name)
		for domain in domains:
			self.predictDomain(domain)
			if None != eagerUpload_conn:
				cursor = eagerUpload_conn.cursor(MySQLdb.cursors.DictCursor)
				self.uploadPredictions(cursor, domain, eager_table_name)
				domain.predictions=Predictions(uploaded=1)

	def uploadBenchmarkPredictions(self, domains, conn, table_name):
		cursor = conn.cursor(MySQLdb.cursors.DictCursor)
		self.createResultsTable(cursor,table_name)

		for domain in domains:
			self.uploadPredictions(cursor, domain, table_name)
		
		cursor.close()

	def createResultsTable(self, cursor, table_name):
		create_qry = """CREATE TABLE IF NOT EXISTS """+table_name+""" LIKE """+self.__TEMPLATE_TABLE_NAME
		cursor.execute(create_qry)

	def uploadPredictions(self,cursor, domain, table_name):
		qry = "INSERT INTO "+table_name+""" 
		(parent_sequence_key, domain_sequence_key, mf_acc, name, p_llr, l_llr, s_llr, pl_llr, ps_llr, ls_llr, pls_llr, base_llr, type, p_acc, l_acc,sccs) 
		VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
		ON DUPLICATE KEY UPDATE 
		parent_sequence_key=VALUES(parent_sequence_key),
		domain_sequence_key=VALUES(domain_sequence_key),
		mf_acc=VALUES(mf_acc),
		name=VALUES(name),
		p_llr=VALUES(p_llr),
		l_llr=VALUES(l_llr),
		s_llr=VALUES(s_llr),
		pl_llr=VALUES(pl_llr),
		ps_llr=VALUES(ps_llr),
		ls_llr=VALUES(ls_llr),
		pls_llr=VALUES(pls_llr),
		base_llr=VALUES(base_llr),
		type=VALUES(type),
		timestamp=NOW(),
		p_acc=VALUES(p_acc),
		l_acc=VALUES(l_acc),
		sccs=VALUES(sccs)
		"""
		insert_list = self.formatPredictions(domain)
		try:
			#print "Uploading", len(insert_list)
			cursor.executemany(qry, insert_list)
		except:
			print qry % insert_list[0]
			raise
	

