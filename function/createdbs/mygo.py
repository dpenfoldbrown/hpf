
import MySQLdb
from hpf.function.term import Terms,Term

def rows2terms(rows, ecode, full=False):
	go_terms = Terms(ecode = ecode)
	for row in rows:
		if full:
			go_terms.append(Term(a=row["acc"],n=row["name"],t=row["term_type"]))
		else:
			#print row["acc"]
			go_terms.append(Term(a=row["acc"]))

	return go_terms

class Mygo():
	def __init__(self, user, passwd, host, e_codes=None, db_name="mygo"):
		self.user = user
		self.passwd = passwd
		self.host = host
		self.db_name = db_name
		self.conn = None
		self.cursor = None
		self.evidence_codes = e_codes

	def get_conn(self):
		return self.conn

	def connect(self):
		self.conn = MySQLdb.connect(host=self.host, user=self.user, passwd=self.passwd, db=self.db_name)
		if self.conn == None:
			raise Exception
		self.cursor = self.conn.cursor(MySQLdb.cursors.DictCursor)

	def isConnected(self):
		if self.conn == None:
			return False
		else:
			return True

	#kdrew: returns a unique set of go terms that adhere to evidence codes for the parameter sequence id list, full stores more information about term such as type and name
	#kdrew: full returns the full description including name and term type
	def getGoTerms(self, seqs, full=False, ancestors=True):
		#kdrew: a bit of a hack to get the list of evidence codes into format so sql won't complain, adds single quote then joins codes with single quote comma single quote and finishes with single quote
		#kdrew: 'IDA','IMP','IEA'
		codes_fmt = "'"+"','".join(self.evidence_codes)+"'"

		if isinstance(seqs,list):
			#kdrew: format sequence list
			seq_fmt = ','.join(map(str,seqs))
		#kdrew: just one sequence id
		else:
			seq_fmt = str(seqs)

		if ancestors:
			#kdrew: this query returns all terms annotated to a specific seq_id with evidence code
			query = "select distinct t.* from gene_product_seq as gps, association as a, term as t, graph_path as gp, evidence as e where gps.gene_product_id = a.gene_product_id and a.term_id = gp.term2_id and gp.term1_id = t.id and e.association_id = a.id and t.is_obsolete = 0 and a.is_not = 0 and e.code in (%s) and gps.seq_id in (%s)" % (codes_fmt,seq_fmt)
		else:
			query = "select distinct t.* from gene_product_seq as gps, association as a, term as t, evidence as e where gps.gene_product_id = a.gene_product_id and a.term_id = t.id and e.association_id = a.id and t.is_obsolete = 0 and a.is_not = 0 and e.code in (%s) and gps.seq_id in (%s)" % (codes_fmt,seq_fmt)

		#print query

		self.cursor.execute(query)
		rows = self.cursor.fetchall ()
	
		return rows2terms(rows,self.evidence_codes,full)

	def get_all_terms(self, full=False):
 		query = "select distinct t.* from term as t where t.is_obsolete = 0 "
		self.cursor.execute(query)
		rows = self.cursor.fetchall ()
		return rows2terms(rows,self.evidence_codes, full)

