
from hpf.hddb.db import url
from sqlalchemy import create_engine, Table, Column, MetaData
from sqlalchemy.sql import select, and_
from sqlalchemy.orm import sessionmaker

engine = create_engine(url+'pwinters')
metadata = MetaData()
metadata.bind = engine
#Session = sessionmaker()
#Session.configure(bind=engine)

p_experiments = Table('p_experiments', metadata, autoload=True)
p_fpsi = Table('p_fpsi', metadata, autoload=True)
p_ffr = Table('p_ffr', metadata, autoload=True)
p_fmcm = Table('p_fmcm', metadata, autoload=True)
p_fconfident = Table('p_fconfident', metadata, autoload=True)
p_mf_all = Table('p_mf_all', metadata, autoload=True)

def calc_extend():
	seq_funclist = {}
	functions_iea = p_mf_all.select().execute()
	for func in functions_iea:
		try:
			seq_funclist[func.sequence_key].append(func.acc)
		except:
			seq_funclist[func.sequence_key] = list()
			seq_funclist[func.sequence_key].append(func.acc)
			
	domains_extended = {}
	domains_novel = {}
	domains_experiment_extended = {}
	domains_experiment_novel = {}
	exps = p_experiments.select().execute()
	for experiment in exps:
		print experiment
		domains_experiment_novel[experiment.name] = 0
		domains_experiment_extended[experiment.name] = 0
		fdomains = select([p_fconfident.c.parent_sequence_key, p_fconfident.c.domain_sequence_key,p_fconfident.c.size],distinct=True).where(p_fconfident.c.experiment_key == experiment.id)
		print fdomains
		fdomains = fdomains.execute()
		for domain in fdomains:
			print domain.parent_sequence_key, domain.domain_sequence_key
			domain_extend = False
			domain_novel = False
			fpredictions = p_fconfident.select().where(and_(p_fconfident.c.parent_sequence_key == domain.parent_sequence_key, p_fconfident.c.domain_sequence_key == domain.domain_sequence_key)).execute()
			for pred in fpredictions:
				print pred
				
				try:
					if len(seq_funclist[pred.parent_sequence_key]) == 0:
						print "novel"
					elif pred.mf_acc in seq_funclist[pred.parent_sequence_key]:
						print "not novel"
					else:
						print "extended prediction"
						try:
							domain_extend = True
							domains_extended[pred.domain_sequence_key].append(pred.mf_acc)
						except KeyError:
							domains_extended[pred.domain_sequence_key] = list()
							domains_extended[pred.domain_sequence_key].append(pred.mf_acc)
				except KeyError:
					print "novel (keyerror)"
					try:
						domain_novel = True
						domains_novel[pred.domain_sequence_key].append(pred.mf_acc)
					except KeyError:
						domains_novel[pred.domain_sequence_key] = list()
						domains_novel[pred.domain_sequence_key].append(pred.mf_acc)
			
			if domain_extend:
				domains_experiment_extended[experiment.name] += domain.size
			if domain_novel:
				domains_experiment_novel[experiment.name] += domain.size

	print domains_experiment_extended
	print domains_experiment_novel

	print "domains_extended:", len(domains_extended.keys())
	print "domains_novel:", len(domains_novel.keys())

	return domains_experiment_novel, domains_experiment_extended


def main():

	organism = {}
	domains = {}
	psiblast = {}
	fr = {}
	mcm = {}
	novel = {}
	extend = {}


	for experiment in p_experiments.select().execute():
		#print experiment
		organism[experiment[1]] = experiment[1]

		#kdrew: initializing
		domains[experiment[1]] = 0
		psiblast[experiment[1]] = 0
		fr[experiment[1]] = 0
		mcm[experiment[1]] = 0
		novel[experiment[1]] = 0
		extend[experiment[1]] = 0

	for experiment in p_fpsi.select().execute():
		psiblast[experiment[1]] = experiment[3]
		domains[experiment[1]] +=  experiment[3]

	for experiment in p_ffr.select().execute():
		try:
			fr[experiment[1]] = experiment[3]
			domains[experiment[1]] +=  experiment[3]
		except KeyError:
			continue

	for experiment in p_fmcm.select().execute():
		try:
			mcm[experiment[1]] = experiment[3]
			domains[experiment[1]] +=  experiment[3]
		except KeyError:
			continue

	domains_experiment_novel, domains_experiment_extended = calc_extend()

	print "Organism,Domains,PDBBlast,FFAS03,Rosetta,Novel,Extended"
	sorted_keys = organism.keys()
	sorted_keys.sort()
	for key in sorted_keys:
		print str(organism[key]) + "," + 		\
			str(domains[key]) + "," +		\
			str(psiblast[key]) + "," +		\
			str(fr[key]) + "," +		\
			str(mcm[key]) + "," +	\
			str(domains_experiment_novel[key]) + "," +		\
			str(domains_experiment_extended[key]) 




if __name__ == "__main__":
	main()

