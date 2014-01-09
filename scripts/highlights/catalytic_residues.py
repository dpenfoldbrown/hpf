
import os,sys
import MySQLdb
import getopt
import ConfigParser

from Bio import SubsMat
from Bio.SubsMat import MatrixInfo
from Bio.pairwise2 import dictionary_match
from Bio.PDB.PDBParser import PDBParser

from hpf.pdb import mammoth
import hpf.scripts.extract_decoys as ed

FIREDB_TABLE_NAME = "firedb.45_filtered_5aug2009"
#FIREDB_TABLE_NAME = "firedb.45_full_5aug2009"
GAP_PENALTY = -0.5

def usage():
	print "config=(filename of configuration file) ,section=(section of configuration file)" 

def main():
	print "Welcome to Find Catalytic Residues on Decoys"

	config_file = "config.cfg"
	config_section = "Default"

	try:
		print "before getopt"
		opts, args = getopt.getopt(sys.argv[1:], "hc:s:", ["help", "config=", "section="]) 
		print "after getopt"
		print opts
		print args
	except getopt.GetoptError:	
		usage()
		sys.exit(2)
	for opt, arg in opts:
		print opt
		print arg
		if opt in ("-h","--help"):
			usage()
			sys.exit()
		elif opt in ("-c","--config"):
			config_file = arg
		elif opt in ("-s","--section"):
			config_section = arg
		else:
			assert False, "unhandled option"

	cr = CatalyticResidue(config_file, config_section)

	if cr.sequence_key == -1:
		cr.runCatalyticResidue()
	else:
		cr.runCatalyticResidueOne(cr.sequence_key)


class CatalyticResidue():
	
	def __init__(self, config_file=None, section=None):
		config = ConfigParser.RawConfigParser()
		config.read(config_file)

		if None == config_file:
			self._default_init()
			return

		#kdrew: blast parameters
		self.hpf_db = config.get(section, 'hpf_db')
		self.hpf_host = config.get(section, 'hpf_host')
		self.hpf_port = config.getint(section, 'hpf_port')


		self.sql_user= config.get(section, 'sql_user')
		self.sql_password= config.get(section, 'sql_password')

		self.tmp_dir= config.get(section, 'tmp_dir')
		#self.decoy_extract_script = config.get(section, 'decoy_extract_script')

		self.high_conf_limit = config.getfloat(section,'high_conf_limit')

		self.exclude_experiments_str = config.get(section,'exclude_experiments')
		self.exclude_experiments = map(int,self.exclude_experiments_str.split(','))

		self.sequence_key = config.getint(section,'sequence_key')

		self.score_threshold = config.getint(section,'score_threshold')
		self.individual_score_threshold = config.getint(section,'individual_score_threshold')


		self.by_random = config.getboolean(section, 'by_random')
		self.rand_seed = config.getint(section, 'rand_seed')
		self.sample_size = config.getint(section, 'sample_size')

		self.submat = dictionary_match(SubsMat.SeqMat(MatrixInfo.blosum62)) 

		self.hpf_conn = MySQLdb.connect(host=self.hpf_host, user=self.sql_user, passwd=self.sql_password, db=self.hpf_db, port=self.hpf_port)

		self.pdbparser = PDBParser(PERMISSIVE=1)

	def runCatalyticResidueOne(self, sequence_key):
		cursor = self.hpf_conn.cursor(MySQLdb.cursors.DictCursor)
		query = """
				SELECT md.sequence_key, md.experiment_astral_ac as astral_id, a.pdb_id, a.chain_id, md.probability, md.cluster_center_index, md.structure_key
				FROM hpf.mcmData AS md, pdb.astral95_1_75 AS a 
				WHERE a.astral_id = md.experiment_astral_ac AND md.sequence_key = %s
			"""

		print query
	
		cursor.execute(query,(sequence_key,))
		high_conf_row = cursor.fetchone()
		self.mapStructureOne(high_conf_row)

	def runCatalyticResidue(self):
		cursor = self.hpf_conn.cursor(MySQLdb.cursors.DictCursor)

		#kdrew: difference from above query is probability greater than instead of sequence_key equal to
		query = """
				SELECT md.sequence_key, md.experiment_astral_ac as astral_id, a.pdb_id, md.probability, md.cluster_center_index, md.structure_key
				FROM hpf.mcmData AS md, pdb.astral95_1_75 AS a 
				WHERE a.astral_id = md.experiment_astral_ac AND md.probability > %s
			"""

		print query
	
		cursor.execute(query,(self.high_conf_limit,))
		high_conf_rows = cursor.fetchall ()

		for high_conf_row in high_conf_rows:
			#high_conf_row = cursor.fetchone()
			self.mapStructureOne(high_conf_row)

	def mapStructureOne(self, high_conf_row):
		#print high_conf_row
		exp_key, exp_name = self.getExperiment(high_conf_row['sequence_key'])
		if exp_key in self.exclude_experiments:
			return
		print exp_name, high_conf_row["sequence_key"], high_conf_row["astral_id"], high_conf_row["pdb_id"], high_conf_row["probability"]
		chain_id = high_conf_row["chain_id"][0]
		pdb_filename = self.downloadPDB(high_conf_row["pdb_id"], tmp_dir = self.tmp_dir)

		decoy_filename = self.extractDecoy(high_conf_row["sequence_key"], self.tmp_dir).values()[0]

		print pdb_filename, decoy_filename
		if not os.path.exists(pdb_filename) or not os.path.exists(decoy_filename): 
			print "error with files existing"
			return

		# Perform Mammoth alignment
		cl = mammoth.MammothCL(experiment=pdb_filename, prediction=decoy_filename)
		mm = mammoth.do_alignment(cl)

		pdb_residues = self.get_pdb_residues(pdb_filename, high_conf_row["pdb_id"], chain_id)
		exp2pred, pdb2exp, pdb2mammoth = self.index_mammoth_alignment(mm, pdb_residues)

		cat_res = self.mapCatalyticResidues(high_conf_row["pdb_id"], high_conf_row["astral_id"][5])
		print "Catalytic Residues"
		self.print_map(mm, cat_res, exp2pred, pdb2exp,pdb2mammoth)

		firedb_res = self.mapFireDbResidues(high_conf_row["pdb_id"], high_conf_row["astral_id"][5])
		print "FireDB Residues"
		self.print_map(mm, firedb_res, exp2pred, pdb2exp,pdb2mammoth)
		

	def get_pdb_residues(self, pdb_filename, pdb_id, chain_id):
		structure = self.pdbparser.get_structure(pdb_id, pdb_filename)
		chain = structure[0][chain_id]
		#print chain.get_id(), " ", len(chain)
		#return chain.get_list()[0].get_id()[1]
		return chain.get_list()

	def print_map(self,mm, important_res, exp2pred, pdb2exp, pdb2mammoth):
		score_total = 0
		conserved_res = []
		identical_res = []
		for res in important_res:
			if len(mm.pred_seq) > pdb2exp[res]:
				if mm.pred_seq[pdb2mammoth[res]] != '.':
					subscore = self.submat(mm.exp_seq[pdb2mammoth[res]], mm.pred_seq[pdb2mammoth[res]])
				else:
					subscore = GAP_PENALTY

				if subscore > self.individual_score_threshold:
					conserved_res.append(res)

				score_total += subscore

				if mm.pred_seq[pdb2mammoth[res]] == mm.exp_seq[pdb2mammoth[res]]:
					print "match: ",res,"->",pdb2mammoth[res],"(",exp2pred[res],") : ", mm.exp_seq[pdb2mammoth[res]] ," : ", mm.pred_seq[pdb2mammoth[res]], " = ",subscore
					identical_res.append(res)
				else:
					print "no match: ",res,"->",pdb2mammoth[res],"(",exp2pred[res],") : ", mm.exp_seq[pdb2mammoth[res]] ," : ", mm.pred_seq[pdb2mammoth[res]], " = ",subscore
				

		pymol_res_select_pdb = "select resi %s" % '+'.join(str(i) for i in important_res)
		print "Ligand binding on exp: ", pymol_res_select_pdb

		pymol_pred_select_pdb = "select resi %s" % '+'.join(str(exp2pred[pdb2exp[i]]) for i in important_res)
		print "Ligand binding on pred: ", pymol_pred_select_pdb

		pymol_pred_conserved_select_pdb = "select resi %s" % '+'.join(str(exp2pred[pdb2exp[i]]) for i in conserved_res)
		print "Ligand binding conserved on pred: ", pymol_pred_conserved_select_pdb

		pymol_pred_identical_select_pdb = "select resi %s" % '+'.join(str(exp2pred[pdb2exp[i]]) for i in identical_res)
		print "Ligand binding identical on pred: ", pymol_pred_identical_select_pdb
		
		print "score_total: ", score_total
		if score_total > self.score_threshold:
			print "**** QUALITY MATCH ****"

		print "\n"


	def mapCatalyticResidues(self, pdb_id, chain_id):
		cat_rows = self.getCatalyticResidues(pdb_id, chain_id)

		#kdrew: map catalytic residues onto mammoth alignment
		#return self.index_mammoth_alignment(mammoth_alignment, [i['residue_number'] for i in cat_res], pdb_residues)
		cat_res = [i['residue_number'] for i in cat_rows]
		return cat_res

	def mapFireDbResidues(self, pdb_id, chain_id):
		firedb_rows = self.getFireDbResidues(pdb_id, chain_id)

		firedb_res = []
		for firedb_row in firedb_rows:
			#print firedb_row
			for i in firedb_row['residue_numbers'].split(','):
				try:
					firedb_res.append(int(i))
				except:
					continue

		#print firedb_res
		#kdrew: map catalytic residues onto mammoth alignment
		#return self.index_mammoth_alignment(mammoth_alignment, firedb_res, pdb_residues)
		return firedb_res

#	def index_mammoth_alignment(self, alignment, cat_res_list, pdb_residues):
#		map_res = {}
#		exp2pred = {}
#		for res in cat_res_list:
#			#kdrew: i is mammoth numbering, j is experiment (pdb) numbering, k is prediction (decoy) numbering
#			j = k = 0
#			
#			for i in xrange(len(alignment.exp_seq)):
#				if alignment.exp_seq[i] != '.':
#					j +=1
#				if alignment.pred_seq[i] != '.':
#					k +=1
#				if res-1 == j:
#					map_res[res] = i+1
#				exp2pred[j+1] = k+1
#				print pdb_residues[j+1].get_resname(), res, i+1, j, pdb_residues[j+1].get_id()[1] ,  k+1
#
#		return map_res, exp2pred

	def index_mammoth_alignment(self, alignment, pdb_residues):
		exp2pred = {}
		pdb2exp = {}
		pdb2mammoth = {}

		#kdrew: i is mammoth numbering, j is experiment (pdb) numbering, k is prediction (decoy) numbering
		j = k = 0

		print pdb_residues[0].get_resname()
		
		for i in xrange(len(alignment.exp_seq)):
			if alignment.exp_seq[i] != '.':
				j +=1
			if alignment.pred_seq[i] != '.':
				k +=1
			exp2pred[j+1] = k+1
			pdb2exp[pdb_residues[j].get_id()[1]] = j+1
			pdb2mammoth[pdb_residues[j].get_id()[1]] = i+1
			print pdb_residues[j].get_resname(), i, j, pdb_residues[j].get_id()[1] ,  k

		return exp2pred, pdb2exp, pdb2mammoth

	def getExperiment(self, sequence_key):
		cursor = self.hpf_conn.cursor(MySQLdb.cursors.DictCursor)
		query = """
				SELECT e.id, e.name 
				FROM hpf.experiment AS e, hpf.domain AS d, hpf.protein AS p
				WHERE e.id = p.experiment_key AND p.sequence_key = d.parent_sequence_key
				AND d.domain_sequence_key = %s
			"""
		cursor.execute(query,(sequence_key,))

		one = cursor.fetchone()
		return (one['id'], one['name'])
	
	def getFireDbResidues(self, pdb_id, chain_id):
		cursor = self.hpf_conn.cursor(MySQLdb.cursors.DictCursor)

		query = """
				SELECT * 
				FROM """+ FIREDB_TABLE_NAME+"""
				WHERE pdb_chain_id = %s 
			"""
		cursor.execute(query,(pdb_id+chain_id,))
		return cursor.fetchall()


	def getCatalyticResidues(self, pdb_id, chain_id):
		cursor = self.hpf_conn.cursor(MySQLdb.cursors.DictCursor)

		query = """
				SELECT * 
				FROM catalytic_site_atlas.CSA_2_2_11 
				WHERE pdb_id = %s 
				AND chain_id = %s
			"""
		cursor.execute(query,(pdb_id,chain_id))
		return cursor.fetchall()

	def extractDecoy(self, sequence_key, tmp_dir='~/tmp/'): 
		ed.opts['silent'] = True
		mcms = ed.mcm([(sequence_key,"")], 1)
		return ed.decoys(mcms, tmp_dir)

	def downloadPDB(self, pdb_id, tmp_dir='~/tmp/', force_download = False):
		filename = tmp_dir + pdb_id + ".pdb"
		filename_gz = filename+".gz"
		log_file = tmp_dir + "wget.log"

		if os.path.exists(filename) and not force_download: 
			return filename

		#kdrew: construct url
		#kdrew: example url ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/02/pdb402d.ent.gz
		base_url = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/%s/pdb%s.ent.gz"

		#kdrew: find middle two chars in pdb_id
		middle2 = pdb_id[1]+pdb_id[2]
		pdb_url = base_url % (middle2, pdb_id)

		wget_cmd = "wget -o %s -O %s %s" % (log_file, filename_gz, pdb_url)
		gunzip_cmd = "gunzip --force %s" % (filename_gz)

		#kdrew: download
		os.system(wget_cmd)
		#kdrew: unzip
		os.system(gunzip_cmd)

		return filename

if __name__ == "__main__":
	main()


