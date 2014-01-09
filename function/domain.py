from structure import *
from term import *
from prediction import *

class Domain:
	def __init__(self,psk=0,dsk=0, conn=None, evidence_codes=DEFAULT_CODES, scale=1.0):
		self.parent_sequence_key = psk
		self.domain_sequence_key = dsk
		self.structures = Structures()
		self.proc_terms = Terms(evidence_codes)
		self.loc_terms = Terms(evidence_codes)
		self.func_terms = Terms(evidence_codes)
		self.predictions = Predictions()
		self.scale = scale
		self.conn = conn
		#kdrew: if a connection is passed in, load in process, localization and function terms
		if(None != conn):
			self.load_process(conn)
			self.load_localization(conn)
			self.load_function(conn)

	def __repr__(self):
		return `self.parent_sequence_key`+" "+`self.domain_sequence_key`+" "+`len(self.structures)`+" "+`len(self.proc_terms)`+" "+`len(self.loc_terms)`+" "+`len(self.func_terms)`+" "+`len(self.predictions)`
	
	def append_prediction(self, p):
		self.predictions.append(p)
	def set_predictions(self, P):
		self.predictions = P

	def get_parent_sequence_key(self):
		return self.parent_sequence_key
	def get_domain_sequence_key(self):
		return self.domain_sequence_key

	def load_process(self, conn,term_table=GO_TABLE,mygo_database=""):
		self.proc_terms.load_terms(conn, self.parent_sequence_key,BIOLOCAL_PROCESS,term_table=term_table,mygo_database=mygo_database)

	def load_localization(self,conn,term_table=GO_TABLE,mygo_database=""):
		self.loc_terms.load_terms(conn, self.parent_sequence_key,CELLULAR_COMPONENT,term_table=term_table,mygo_database=mygo_database)

	def load_function(self,conn,term_table=GO_TABLE,mygo_database=""):
		self.func_terms.load_terms(conn, self.parent_sequence_key,MOLECULAR_FUNCTION,term_table=term_table,mygo_database=mygo_database)
		
	def load_terms(self,conn,term_table=GO_TABLE,mygo_database=""):
		self.load_localization(conn,term_table=term_table,mygo_database=mygo_database)
		self.load_process(conn,term_table=term_table,mygo_database=mygo_database)
		self.load_function(conn,term_table=term_table,mygo_database=mygo_database)
		
	def scale_structures(self):
		self.structures.scale_structures(self.scale)


class Domains(list):
	"holds list of domains"

	def __init__(self, conn, data_table, evidence_codes, scale=1.0):
		#self.benchmark = benchmark
		self.conn = conn
		self.scale = scale
		self.data_table = data_table
		self.evidence_codes = evidence_codes

	def load_domains(self, typeOfDomain):
		#kdrew: import domains from database
		cursor = self.conn.cursor(MySQLdb.cursors.DictCursor)
		qry = "SELECT DISTINCT parent_sequence_key, domain_sequence_key FROM "+self.data_table
		print(qry);
		cursor.execute(qry)
	
		rows = cursor.fetchall ()
		for row in rows:
			d_tmp = apply(typeOfDomain, (row["parent_sequence_key"],row["domain_sequence_key"], self.conn, self.evidence_codes, self.scale))
			#d_tmp.load_mcm_structures(self.conn,self.scale, self.data_table)
			self.append(d_tmp)

class FRDomain(Domain):
	def __init__(self,psk=0,dsk=0, conn=None,evidence_codes=DEFAULT_CODES, scale=0.9):
		Domain.__init__(self,psk,dsk,conn,evidence_codes,scale)
		if(None != conn):
			self.load_structures(conn, scale)
		self.type = "fr"
	
	def load_structures(self,conn, scale=0.9):
		self.structures = FR_Structures()
		self.structures.load_structures(conn, self.domain_sequence_key)
		self.structures.combine_multiple_sf()
		self.structures.scale_structures(scale=scale)

class FRDomains(Domains):
	"holds list of fr domains"

	def __init__(self, conn,  data_table="fr_data", evidence_codes=DEFAULT_CODES, scale=0.9):
		Domains.__init__(self, conn, data_table, evidence_codes, scale)
		self.load_domains()
		
	def load_domains(self):
		Domains.load_domains(self,fr_domain.FRDomain)
		


class MCMDomain(Domain):
	def __init__(self,psk=0,dsk=0, conn=None, evidence_codes=DEFAULT_CODES, scale=0.5):
		Domain.__init__(self,psk,dsk,conn,evidence_codes,scale)
		if(None != conn):
			self.load_structures(conn, scale)
		self.type = "mcm"
	
	def load_structures(self,conn, scale=1.0):
		self.structures = MCM_Structures()
		self.structures.load_structures(conn, self.domain_sequence_key)
		self.structures.combine_multiple_sf()
		self.structures.scale_structures(scale=scale)

	def store_mcm_structures(self,conn, drop=1, data_table="mcm_data_adjusted"):
		self.structures.store_mcm_structures(conn, self.domain_sequence_key, self.parent_sequence_key, drop, data_table)

class MCMDomains(Domains):
	"holds list of mcm domains"

	def __init__(self, conn,  data_table="mcm_data", evidence_codes=DEFAULT_CODES, scale=0.5):
		Domains.__init__(self, conn, data_table, evidence_codes, scale)
		self.load_domains()
		
	def load_domains(self):
		Domains.load_domains(self,mcm_domain.MCMDomain)
	
class PSIDomain(Domain):
	def __init__(self,psk=0,dsk=0, conn=None, evidence_codes=DEFAULT_CODES, scale=1.0):
		Domain.__init__(self,psk,dsk,conn,evidence_codes,scale)
		if(None != conn):
			self.load_structures(conn, scale)
		self.type = "psi"

	def load_structures(self,conn, scale=1.0):
		self.structures = PSI_Structures()
		self.structures.load_structures(conn, self.domain_sequence_key)
		self.structures.combine_multiple_sf()
		self.structures.scale_structures(scale=scale)
		
class PSIDomains(Domains):
	"holds list of psi domains"

	def __init__(self, conn,  data_table="psi_data", evidence_codes=DEFAULT_CODES, scale=1.0):
		Domains.__init__(self, conn, data_table, evidence_codes, scale)
		self.load_domains()
		
	def load_domains(self):
		Domains.load_domains(self,psi_domain.PSIDomain)

class NewAnnotatedDomain(Domain):
	def __init__(self,psk=0,dsk=0, evidence_codes=DEFAULT_CODES):
		self.parent_sequence_key = psk
		self.domain_sequence_key = dsk
		self.evidence_codes = evidence_codes
		self.old_functions = Terms()
		self.new_functions = Terms()

	def set_old_functions(self,old_terms):
		self.old_functions = old_terms
	def get_old_functions(self):
		return self.old_functions

	def set_new_functions(self,new_terms):
		self.new_functions = new_terms
	def get_new_functions(self):
		return self.new_functions

class NewAnnotatedDomains(Domains):
	"holds list of newly annotated domains"

	def __init__(self, conn, seq_table, domain_table, evidence_codes):
		self.conn = conn
		self.seq_table = seq_table
		self.domain_table = domain_table
		self.evidence_codes = evidence_codes

	def load_domains(self):
		#kdrew: import domains from database
		cursor = self.conn.cursor(MySQLdb.cursors.DictCursor)
		qry = "SELECT DISTINCT d.parent_sequence_key, d.domain_sequence_key FROM "+self.seq_table+" AS s, "+self.domain_table+" AS d WHERE d.parent_sequence_key = s.sequence_key"
		print(qry);
		cursor.execute(qry)
	
		rows = cursor.fetchall ()
		for row in rows:
			d_tmp = NewAnnotatedDomain(row["parent_sequence_key"],row["domain_sequence_key"], self.evidence_codes)
			self.append(d_tmp)


