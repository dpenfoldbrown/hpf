
import MySQLdb
from hpf.function.predictor import Predictor

class AstralList(list):
	def __init__(self, a_db):
		self.astral_db = a_db

class AstralEntry(Predictor):
	def __init__(self, scop_id, sccs_id, header=None, sequence=None):
		self.scop_id = scop_id
		self.sccs_id = sccs_id
		self.header = header
		self.sequence = sequence

	def get_id(self):
		return self.sccs_id

class Astral():
	def __init__(self, user, passwd, host, db_name="hpf", table="astral"):
		self.user = user
		self.passwd = passwd
		self.host = host
		self.db_name = db_name
		self.table = table
		self.conn = None
		self.cursor = None

	def connect(self):
		self.conn = MySQLdb.connect(host=self.host, user=self.user, passwd=self.passwd, db=self.db_name)
		self.cursor = self.conn.cursor(MySQLdb.cursors.DictCursor)

	def isConnected(self):
		if self.conn == None:
			return False
		else:
			return True

	def get_all_superfamilies(self):
 		query = "select distinct SUBSTRING_INDEX(sccs,'.',3) as superfamily from "+self.table+" as a"
		self.cursor.execute(query)
		rows = self.cursor.fetchall ()
		superfamilies = []

		for row in rows:
			superfamilies.append(row["superfamily"])

		return superfamilies

class SuperfamilyEntry(AstralEntry):
	def __init__(self, scop_id, sccs_id, header=None, sequence=None):
		AstralEntry.__init__(self,scop_id,sccs_id,header,sequence)
		self.superfamily = self.sccs_id.rpartition('.')[0]
	def __init__(self,sf_id):
		self.superfamily = sf_id
	def get_id(self):
		return self.superfamily




