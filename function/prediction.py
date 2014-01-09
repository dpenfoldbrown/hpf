
import MySQLdb

from term import MOLECULAR_FUNCTION, BIOLOCAL_PROCESS, CELLULAR_COMPONENT, Term

class Prediction:
	def __init__(self,t=None, proc_feature = None, loc_feature = None,
		s_score = 0, l_score = 0, p_score = 0,
		pl_score = 0, ps_score = 0, ls_score = 0,
		pls_score = 0, base_score = 0,
		st=None):

		self.function = t
		self.proc_feature = proc_feature
		self.loc_feature =loc_feature 
		self.s_score = s_score
		self.l_score = l_score
		self.p_score = p_score
		self.pl_score = pl_score
		self.ps_score = ps_score
		self.ls_score = ls_score
		self.pls_score = pls_score
		self.base_score = base_score
		self.score_type = st

	#kdrew: implement __eq__
	def __eq__(self, other):
		return 0 == cmp(self.function.get_id(),other.get_id())

	def __repr__(self):
		return `self.function` +" "+ `self.proc_feature` +" "+ `self.loc_feature` +"\n"+ "s:%f l:%f p:%f pl:%f ps:%f ls:%f pls:%f base:%f" % (self.s_score, self.l_score, self.p_score, self.pl_score, self.ps_score, self.pl_score, self.pls_score, self.base_score) +" "+ `self.score_type`

	def update_proc_feature(self,p):
		self.proc_feature = p

	def update_loc_feature(self,l):
		self.loc_feature = l

	def update_s_score(self,s):
		self.s_score = s
	def update_l_score(self,l):
		self.l_score = l
	def update_p_score(self,p):
		self.p_score = p
	def update_ps_score(self,ps):
		self.ps_score = ps
	def update_ls_score(self,ls):
		self.ls_score = ls
	def update_pl_score(self,pl):
		self.pl_score = pl
	def update_pls_score(self,pls):
		self.pls_score = pls
	def update_base_score(self,base):
		self.base_score= base

class Predictions(list):
	"list of predictions"

	def __init__(self, uploaded=0, conn=None, table=None):
		self.uploaded = uploaded
		self.conn = conn
		self.table_name = table

	def __repr__(self):
		ret_string = ""
		for pred in self:
			ret_string += pred.__repr__()
			ret_string += "\n"
		return ret_string

	#kdrew: returns a prediction when given a go acc
	def getByAcc(self,acc):
		try:
			return self[self.index(Term(acc,None,None))]
		except ValueError:
			return False

	#kdrew: filter predictions by minimum llr score, base_llr is exception (filters by maximum)
	def filter(self, s_llr=-99, l_llr=-99, p_llr=-99, ps_llr=-99, pl_llr=-99, ls_llr=-99, pls_llr=-99, base_llr=99):
		ret_preds = Predictions()
		for pred in self:
			if (	 pred.s_score >= s_llr 
					and pred.l_score >= l_llr 
					and pred.p_score >= p_llr 
					and pred.ps_score >= ps_llr 
					and pred.ls_score >= ls_llr 
					and pred.pl_score >= pl_llr 
					and pred.pls_score >= pls_llr 
					and pred.base_score <= base_llr
				):

				ret_preds.append(pred)

		return ret_preds

	def load_predictions(self, psk, dsk):
		cursor = self.conn.cursor(MySQLdb.cursors.DictCursor)

		query = """SELECT * FROM """+self.table_name+""" AS pds 
				WHERE pds.parent_sequence_key = %s AND pds.domain_sequence_key = %s
			"""

		cursor.execute(query, (psk, dsk))
		pred_rows = cursor.fetchall()
		for pred_row in pred_rows:
			mf_term = Term(pred_row["mf_acc"], pred_row["name"], MOLECULAR_FUNCTION)
			bp_term = Term(pred_row["p_acc"], BIOLOCAL_PROCESS)
			loc_term = Term(pred_row["l_acc"], CELLULAR_COMPONENT)
			pred = Prediction(t=mf_term, proc_feature = bp_term, loc_feature = loc_term,
						s_score = pred_row["s_llr"], l_score = pred_row["l_llr"], p_score = pred_row["p_llr"],
						pl_score = pred_row["pl_llr"], ps_score = pred_row["ps_llr"], ls_score = pred_row["ls_llr"],
						pls_score = pred_row["pls_llr"], base_score = pred_row["base_llr"],
						st=pred_row["type"])

			self.append(pred)


