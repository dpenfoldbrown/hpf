
from hpf.hddb.db import url
from sqlalchemy import create_engine, Table, Column, MetaData
from sqlalchemy.orm import sessionmaker

engine = create_engine(url+'since_solved')
metadata = MetaData()
metadata.bind = engine

p_table = Table('p_sap_table', metadata, autoload=True)

def main():

	#select * from p_mcm_redux as pmr, hpf.sequence as s where s.id = pmr.sequence_key order by SUBSTRING_INDEX(pred_sf,'.',1),top_prob DESC

	since_solved = p_table.select().execute()
	print "sequence_key,predicted superfamily,true superfamily,mcm score,correct,sequence"
	for entry in since_solved:
		print str(entry["sequence_key"]) + "," +	\
			str(entry["pred_sf"]) + "," +		\
			str(entry["true_sf"]) + "," +		\
			str(entry["mcm_score"]) + "," +		\
			str(entry["correct"]) + "," +		\
			str(entry["sequence"])
			



if __name__ == "__main__":
	main()
