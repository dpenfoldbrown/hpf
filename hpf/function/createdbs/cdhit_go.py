
import os
import hpf.cdhit
from hpf.function.predictor import Predictor

class CDHitGO():
	
	def __init__(self, ch_e, i_f, o_f=None, id_c=.9, len_c=0):
		self.identity_cutoff = id_c
		self.length_cutoff = len_c
		self.in_file = i_f
		self.cd_hit_exe = ch_e
		if o_f != None:
			self.out_file = o_f
		else:
			self.out_file = self.in_file.rpartition('.')[0]+".cd_hit.out"

		self.clstr_file = self.out_file+".clstr"

	def runCDHit(self): 
		cmd = self.cd_hit_exe+" -i "+self.in_file+" -o "+self.out_file+" -c "+str(self.identity_cutoff)+" -s "+str(self.length_cutoff)
		print cmd
		os.system(cmd)

	def printCDHit(self):
		parser = GO_CDHitParser(self.clstr_file)
		for cluster in parser:
			print "Cluster",cluster.id
			for entry in cluster:
				print "\t",entry.id,entry.perc,entry.repr


	#kdrew: now deprecated
	def getClusterGoTerms(self, mygo):
		clstr_dict = {}
		parser = GO_CDHitParser(self.clstr_file)
		for cluster in parser:
			clstr_terms = []
			seq_ids = []
			clstr_rep = 0
			for entry in cluster:
				seq_ids.append(entry.id)
				if entry.repr:
					clstr_rep = entry.id
		
			clstr_terms = mygo.getGoTerms(seq_ids)
			print "id: ",clstr_rep, " num of terms: ", len(clstr_terms)
			clstr_dict[clstr_rep] = clstr_terms

		return clstr_dict

	#kdrew: add terms and superfamiles to a list of predictors for a given cluster representative
	def getClusters(self,filtered_seqs, mygo, freq_metric=None, mustHaveGO = True):
		clstr_dict = {}
		parser = GO_CDHitParser(self.clstr_file)
		for cluster in parser:
			clstr_accs = []
			seq_ids = []
			clstr_rep = 0
			hasRemove = False
			for entry in cluster:
				if entry.id == "remove":
					hasRemove = True
				else:
					seq_ids.append(entry.id)
				if entry.repr:
					clstr_rep = entry.id

			#kdrew: cluster must not have "remove" sequence
			if hasRemove:
				print "cluster has remove: continue"
				continue
		
			clstr_accs = mygo.getGoTerms(seq_ids)
			#kdrew: if there are no go terms then go to the next iteration, must have atleast one go term to be included
			if mustHaveGO and len(clstr_accs) == 0:
				continue

			clstr_accs.append(Predictor(filtered_seqs[clstr_rep].superfamily))
			print "id: ",clstr_rep, " num of terms: ", len(clstr_accs)

			if freq_metric == None:
				clstr_dict[clstr_rep] = clstr_accs
			else:
				freq_metric.compute_metric_cluster(clstr_accs)
			

		return clstr_dict


class GO_CDHitParser(hpf.cdhit.CDHitParser):

	def next(self):
		#>Cluster 0
		#0       459aa, >3|Q9XHP0|11S glob... *

		if(self._finished()):
			raise StopIteration()
		line = self._get()
		if not line.startswith(">"):
			raise Exception("Cluster",line)
		cluster = hpf.cdhit.CDHitCluster(int(line.split()[-1]))
		while True:
			if self._peek().startswith(">") or self._finished():
				break
			line = self._get().replace("\t"," ")
			num = line.split()[0]
			aa = line.split()[1].split('aa')[0]
			id = line.split('>')[1].split('|')[0]
			if line.find('%') > -1:
				perc = int(line.split('... at')[1].strip().strip('%'))
				repr = False
			elif line.find('*') > -1:
				perc = 100
				repr = True
			else:
				raise Exception("Entry",line)
			cluster.append(hpf.cdhit.CDHitEntry(id,aa,num,perc,repr))
		return cluster



