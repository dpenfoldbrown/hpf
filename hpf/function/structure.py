from predictor import *

class Structure(Predictor):
	
	def __init__(self,sf, p,domain_type=None):
		self.superfamily = str(sf)
		self.probability = float(p)
		self.domain_type = domain_type

	#kdrew: predictor interface
	def get_id(self):
		return self.superfamily

	def duplicate(self):
		return Structure(self.superfamily, self.probability)

	def update_prob(self, p):
		self.probability = p

	def get_prob(self):
		return self.probability

	def __repr__(self):
		return `self.superfamily`+" "+`self.probability`
	
	def __cmp__(self, other):
		return cmp(self.probability, other.probability)
	
	def __eq__(self, other):
		return self.probability==other.probability and self.superfamily==other.superfamily

def compare(struct1, struct2):
	if struct1.probability < struct2.probability:
		return 1
	elif struct1.probability == struct2.probability:
		return 0
	else: #struct1.probability > struct2.probability:
		return -1

class Structures(list):

	#1-(1-probability)*(1-probability2)
	#kdrew: combine probabilities by the above formula
	#remove structures as we come to them, find 
	def combine_multiple_sf(self):
		# Create a copy
		structs = []
		structs.extend(self)
		#for struct in structs[:]:
		while 0 < len(structs):
			struct = structs.pop(0)
			#print "struct: ",struct
			multiply_prob = (1-struct.get_prob())

			for struct2 in structs[:]:
				#print "struct2: ",struct2
				if struct.superfamily == struct2.superfamily:
					multiply_prob = multiply_prob*(1-struct2.get_prob())
					structs.pop(structs.index(struct2))
					self.pop(self.index(struct2))

			self[self.index(struct)].update_prob(1-multiply_prob)
			#print "updated_struct: ",self[self.index(struct)]

	def scale_structures(self, scale=1.0):
		#print self
		for struct in self:
			struct.update_prob(scale*(struct.get_prob()))
		
#		sum_prob = float(sum([float(struct.get_prob()) for struct in self]))
#		for struct in self:
#			struct.update_prob(float(scale*(struct.get_prob()/sum_prob)))
		#struct.update_prob(scale*(struct.get_prob()))
		#print self
				
	def clear_list(self):
		while 0 < len(self):
			self.pop(0)
#
class PSI_Structure(Structure):
	"holds one psiblast structure"	
	def __init__(self,dk,sf,p):
		Structure.__init__(self,sf,p)
		self.domainKey = dk

	def duplicate(self):
		return PSI_Structure(self.domainKey, self.superfamily, self.probability)

class PSI_Structures(Structures):
	"holds multiple psiblast structure predictions"

	def load_structures(self, conn, dsk, clear=1, data_table="psi_data"):
		if clear:
			self.clear_list()

		#kdrew: import structures from database
		cursor = conn.cursor(MySQLdb.cursors.DictCursor)
		qry = "select * from "+data_table+" where domain_sequence_key = "+str(dsk)
		cursor.execute(qry)
	
		rows = cursor.fetchall ()
		for row in rows:
			self.append(psi_structure.PSI_Structure(row["domain_key"],row["sccs"],float(row["probability"])))

		cursor.close()

	
class MCM_Structure(Structure):
	"holds one mcm structure"	
	def __init__(self,mdk,sf,p):
		Structure.__init__(self,sf,p)
		self.mcmDataKey = mdk

	def duplicate(self):
		return MCM_Structure(self.mcmDataKey, self.superfamily, self.probability)
	
		#kdrew: implement __eq__
	def __eq__(self, other):
		return 0 == cmp(self.mcmDataKey,other.mcmDataKey)

class MCM_Structures(Structures):
	"holds multiple mcm structure predictions"

	def load_structures(self, conn, dsk, psk=None, clear=1, data_table="mcm_data"):
		self.load_mcm_structures(conn=conn, dsk=dsk, psk=psk, clear=clear, data_table=data_table)

	def load_mcm_structures(self, conn, dsk, psk=None, clear=1, data_table="mcm_data"):
		if clear:
			self.clear_list()

		#kdrew: import structures from database
		cursor = conn.cursor(MySQLdb.cursors.DictCursor)
		qry = "select * from "+data_table+" where domain_sequence_key = "+str(dsk)
		if psk != None:
			qry = qry+" and parent_sequence_key = "+str(psk)
		cursor.execute(qry)
	
		rows = cursor.fetchall ()
		for row in rows:
			#self.structures.append(Structure(row[0],row[3],row[4]))
			self.append(mcm_structure.MCM_Structure(row["mcmdata_key"],row["sccs"],float(row["probability"])))

		cursor.close()

	def store_mcm_structures(self, conn, dsk, psk, drop=1, data_table="mcm_data_adjusted"):
		cursor = conn.cursor(MySQLdb.cursors.DictCursor)

		if drop:
			drop_qry = "drop table if exists "+data_table
			cursor.execute(drop_qry)

		create_qry = "create table if not exists "+data_table+"( mcmdata_key int(11), parent_sequence_key int(11), domain_sequence_key int(11), sccs varchar(15), probability double, KEY parent_sequence_key (parent_sequence_key), KEY domain_sequence_key (domain_sequence_key))"
		cursor.execute(create_qry)

		insert_qry = "insert into "+data_table+" (mcmdata_key, parent_sequence_key, domain_sequence_key, sccs, probability) VALUES "
		insert_list = []
		for struct in self:
			insert_list.append("(" +str(struct.mcmDataKey)+ "," +str(psk)+ "," +str(dsk)+ ",\'" +struct.superfamily+ "\'," +str(struct.probability)+ ")")
			#print insert_qry+"\n"

		insert_qry = insert_qry+','.join(insert_list)

		cursor.execute(insert_qry)

		cursor.close()

	
class FR_Structure(Structure):
	"holds one fr structure"	
	def __init__(self,dk,sf,p):
		Structure.__init__(self,sf,p)
		self.domainKey = dk

	def duplicate(self):
		return FR_Structure(self.domainKey, self.superfamily, self.probability)

#def compare(struct1, struct2):
#	return cmp(struct1.probability, struct2.probability)



