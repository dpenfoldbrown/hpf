
import MySQLdb
import string
from predictor import *

MOLECULAR_FUNCTION = "molecular_function"
BIOLOCAL_PROCESS = "biological_process"
CELLULAR_COMPONENT = "cellular_component"
#DEFAULT_CODES = ['TAS','IDA','IMP','IGI','IPI','ISS','IEA','IC','IEP','NAS','ND','NR','RCA']
DEFAULT_CODES = ['TAS','IDA','IMP','IGI','IPI','ISS','IEA','hpf_IEA']
GO_TABLE = "hddb_IEA.hddb_IEA_goLite_062009"
class Term(Predictor):
	def __init__(self, a, n=None, t=None, c=None):
		self.acc = a
		self.name = n
		self.type = t
		self.code = c

	#kdrew: predictor interface
	def get_id(self):
		return self.acc

	def get_name(self):
		return self.name

	def __eq__(self,other):
		if hasattr(other, "get_id"):
			return self.get_id() == other.get_id()
		else:
			return False

	def __str__(self):
		return self.get_id()

	def __hash__(self):
		return hash(self.get_id())

	def __repr__(self):
		return `self.acc`+" "+`self.name`+" "+`self.type`+" "+`self.code`

def is_type(go_acc, type, conn):
	query = """select term_type from term where acc = %s"""
	cursor = conn.cursor(MySQLdb.cursors.DictCursor)
	cursor.execute(query,(go_acc,))
	row = cursor.fetchone()
	if None != row:
		return row["term_type"] == type
	else:
		return False
def get_name(go_acc, conn):
	query = """select name from term where acc = %s"""
	cursor = conn.cursor(MySQLdb.cursors.DictCursor)
	cursor.execute(query,(go_acc,))
	row = cursor.fetchone()
	if None != row:
		return row["name"] 
	else:
		return ""

class Terms(list):
	def __init__(self, ecode = DEFAULT_CODES ):
		self.evidence_codes = ecode

	def load_terms(self, conn, psk, go_type, term_table=GO_TABLE, mygo_database=""):
		"""
		Loads all terms and ancestors from the GO DAG by joining known annotations
		with the graph path.
		@param term_table: The table to load sequence annotations from.
		@param mygo_database: The name of the mygo database to use for joining
			terms with the graph_path.  To get implied ancestor terms.
			Defaults to '' which will assume the database is specified by the
			connection.
		"""
		if not mygo_database.endswith(".") and mygo_database != "":
			mygo_database+="."

		#kdrew: allows parameter go_type to be passed in as list or string
		if not isinstance(go_type, list):
			go_type = [go_type,]

		cursor = conn.cursor(MySQLdb.cursors.DictCursor)

		# join up on graph path to get implied ancestral terms.
		query = """
				SELECT DISTINCT t1.acc,t1.name,t1.term_type FROM
				%s g JOIN %sterm t2 ON g.acc=t2.acc AND t2.term_type IN ('%s') AND t2.is_obsolete=0
				JOIN %sgraph_path p ON t2.id=p.term2_id
				JOIN %sterm t1 ON p.term1_id=t1.id AND t1.term_type IN ('%s') AND t1.is_obsolete=0
				WHERE g.evidence_code IN ('%s')
				AND g.sequence_key=%i
				""" % (term_table,mygo_database,string.join(go_type,"','"),mygo_database,mygo_database,string.join(go_type,"','"),string.join(self.evidence_codes,"','"),psk)

		#print query
		cursor.execute(query)
		for row in cursor.fetchall():
			self.append(Term(row["acc"], 
							 row["name"], 
							 row["term_type"]))

	def print_terms(self):
		for term in self:
			#term.print_term()
			print term

	def getFeaturedTerm(self, function, feature_metric):
		return feature_metric.max(function,self)
	
class ProcTerm(Term):
        "holds one biological process term"
	
class LocTerm(Term):
        "holds one cellular component (localization) term"

class FuncTerm(Term):
        "holds one molecular function term"

class FuncTerms(Terms):
	"container for storing all molecular function terms"
	def __init__(self, conn=None,query=None):
		if None != conn:
			self.load_terms(conn, query)

	def load_terms(self, conn,query = None):
		cursor = conn.cursor(MySQLdb.cursors.DictCursor)
		if query == None:
			query = "select DISTINCT t.acc, t.name, t.term_type"
			query = query+ " from term as t"
			query = query+ " where t.term_type = 'molecular_function' and t.is_obsolete=0 "

		print query
		cursor.execute(query)
		rows = cursor.fetchall ()
		for row in rows:
			self.append(FuncTerm(row["acc"], row["name"], row["term_type"]))
		cursor.close()



