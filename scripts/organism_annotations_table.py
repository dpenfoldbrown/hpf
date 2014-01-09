
from hpf.hddb.db import url
from sqlalchemy import create_engine, Table, Column, MetaData
from sqlalchemy.orm import sessionmaker

engine = create_engine(url+'pwinters')
metadata = MetaData()
metadata.bind = engine
#Session = sessionmaker()
#Session.configure(bind=engine)

p_domains = Table('p_domains', metadata, autoload=True)
p_psi = Table('p_psi', metadata, autoload=True)
p_fold = Table('p_fold', metadata, autoload=True)
p_pfam = Table('p_pfam', metadata, autoload=True)
p_msa = Table('p_msa', metadata, autoload=True)
p_unassigned = Table('p_unassigned', metadata, autoload=True)
p_rosetta = Table('p_rosetta', metadata, autoload=True)
p_rosetta_8 = Table('p_rosetta_8', metadata, autoload=True)
p_rosetta_9 = Table('p_rosetta_9', metadata, autoload=True)
p_executable = Table('p_executable', metadata, autoload=True)

def main():

	organism = {}
	proteins = {}
	domains = {}
	psiblast = {}
	fr = {}
	pfam = {}
	msa = {}
	unassigned = {}
	rosetta = {}
	rosetta_8 = {}
	rosetta_9 = {}
	executable = {}

	for experiment in p_domains.select().execute():
		#print experiment
		organism[experiment[1]] = experiment[1]
		proteins[experiment[1]] = experiment[2]
		#domains[experiment[1]] = experiment[3]

		#kdrew: initializing
		domains[experiment[1]] = 0
		psiblast[experiment[1]] = 0
		fr[experiment[1]] = 0
		pfam[experiment[1]] = 0
		msa[experiment[1]] = 0
		unassigned[experiment[1]] = 0
		rosetta[experiment[1]] = 0
		rosetta_8[experiment[1]] = 0
		rosetta_9[experiment[1]] = 0
		executable[experiment[1]] = ''

	for experiment in p_psi.select().execute():
		psiblast[experiment[1]] = experiment[3]
		domains[experiment[1]] +=  experiment[3]

	for experiment in p_fold.select().execute():
		fr[experiment[1]] = experiment[3]
		domains[experiment[1]] +=  experiment[3]

	for experiment in p_pfam.select().execute():
		pfam[experiment[1]] = experiment[3]
		domains[experiment[1]] +=  experiment[3]

	for experiment in p_msa.select().execute():
		msa[experiment[1]] = experiment[3]
		domains[experiment[1]] +=  experiment[3]

	for experiment in p_unassigned.select().execute():
		unassigned[experiment[1]] = experiment[3]
		domains[experiment[1]] +=  experiment[3]

	for experiment in p_rosetta.select().execute():
		rosetta[experiment[1]] = experiment[3]

	for experiment in p_rosetta_8.select().execute():
		rosetta_8[experiment[1]] = experiment[3]

	for experiment in p_rosetta_9.select().execute():
		rosetta_9[experiment[1]] = experiment[3]

	for experiment in p_executable.select().execute():
		executable[experiment[1]] = experiment[2]

	print "Organism,Proteins,Domains,PDBBlast,FFAS03,Pfam,MSA,Heuristic,Rosetta (total),Rosetta (medium conf),Rosetta (high conf),High Resolution Rosetta"
	sorted_keys = organism.keys()
	sorted_keys.sort()
	for key in sorted_keys:
		print str(organism[key]) + "," + 		\
			str(proteins[key]) + "," +		\
			str(domains[key]) + "," +		\
			str(psiblast[key]) + "," +		\
			str(fr[key]) + "," +		\
			str(pfam[key]) + "," +		\
			str(msa[key]) + "," +		\
			str(unassigned[key]) + "," +	\
			str(rosetta[key]) + "," +		\
			str(rosetta_8[key]) + "," +		\
			str(rosetta_9[key]) + "," +		\
			str(executable[key]) 




if __name__ == "__main__":
	main()

