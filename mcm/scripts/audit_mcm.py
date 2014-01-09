from hpf.hddb.db import RosettaConvergence, RosettaCluster, FilesystemOutfile
from hpf.hddb.db import Structure
from hpf.hddb.db import McmData, Session,url

from numpy import average

from sqlalchemy.sql import func, and_
import sqlalchemy
import getopt
import sys

engine = sqlalchemy.create_engine(url+'pwinters')
metadata = sqlalchemy.MetaData(engine)
pdomainsccs = sqlalchemy.Table('p_domain_sccs',metadata, autoload=True)

session = Session()


def usage():
        print "pickle=(sequences to task.pickle , default=None) , pickle_file=(pickle to file, default=tasks.pickle, sequences=(sequences to audit, 'all' or 'paper', default=paper)"


def main():
	sequences = "paper"
	pickle = None
	pickle_file = "tasks.pickle"

        try:
                print "before getopt"
                opts, args = getopt.getopt(sys.argv[1:], "hs:p:f:", ["help", "sequences=", "pickle=","pickle_file="])
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
                elif opt in ("-p","--pickle"):
                        pickle = arg
                elif opt in ("-f","--pickle_file"):
                        pickle_file = arg
                elif opt in ("-s","--sequences"):
                        sequences = arg
                else:
                        assert False, "unhandled option"

	#sequence_keys = [262258,284577]

	if sequences == 'paper':
		sks, nsk = get_paper_sequences()
	elif sequences == 'all':
		sks, nsk = get_all_denovo_sequences()
	else:
		try:
			sks = []
			sks.append(int(sequences))
			nsk = len(sks)
		except:
			print "error in setting sequences"
			os.exit()

	print sks
	audit_dict = run_audit(sks,nsk)
	
	if pickle != None:
		import cPickle
		with open(pickle_file,"w") as handle:
			cPickle.dump(audit_dict[pickle],handle)

def get_all_denovo_sequences():
	sequence_keys = session.query(McmData.sequence_key).distinct()
	num_sequence_keys = session.query(McmData.sequence_key).distinct().count()
	return sum(sequence_keys,()), num_sequence_keys

def get_paper_sequences():
	squery = sqlalchemy.select([pdomainsccs.c.domain_sequence_key],
					sqlalchemy.and_((~(pdomainsccs.c.domain_type.in_(("fold_recognition","psiblast")))), pdomainsccs.c.sccs != "NULL"),
					distinct = True)
	sequence_keys = squery.execute()

	ret_seq_keys = []
	for i in sequence_keys:
		ret_seq_keys.append(i[0])
	
	squery2 = sqlalchemy.select(columns=[sqlalchemy.func.count(pdomainsccs.c.domain_sequence_key)],
					whereclause=sqlalchemy.and_((~(pdomainsccs.c.domain_type.in_(("fold_recognition","psiblast")))),pdomainsccs.c.sccs != "NULL"),
					distinct = True)
	num_sequence_keys = squery2.execute().fetchone()[0]

	return ret_seq_keys, num_sequence_keys

def run_audit(sequence_keys, num_sequence_keys):
	sum_audit_true = 0
	sum_mcm_true = 0
	sum_probability_true = 0
	sum_convergence_true = 0
	sum_clusters_true = 0
	complete_seqs = []
	incomplete_seqs = []
	for sequence_key in sequence_keys:
		print "sequence_key: ", sequence_key
		mcm_flag,probability_flag,convergence_flag,clusters_flag = audit_sequence(sequence_key)
		print "mcm_flag: ", mcm_flag, "probability_flag: ", probability_flag, " convergence_flag: ",convergence_flag, " clusters_flag: ", clusters_flag
		audit_flag = True if (mcm_flag and probability_flag and convergence_flag and clusters_flag) else False
		print "Complete upload: ", audit_flag

		sum_audit_true += 1 if audit_flag else 0
		sum_mcm_true += 1 if mcm_flag else 0
		sum_probability_true += 1 if probability_flag else 0
		sum_convergence_true += 1 if convergence_flag else 0
		sum_clusters_true += 1 if clusters_flag else 0
		
		print "MCM true: ", sum_mcm_true, " / ", num_sequence_keys, "(", 1.0*sum_mcm_true/num_sequence_keys, ")"
		print "Probability true: ", sum_probability_true, " / ", num_sequence_keys, "(", 1.0*sum_probability_true/num_sequence_keys, ")"
		print "Convergence true: ", sum_convergence_true, " / ", num_sequence_keys, "(", 1.0*sum_convergence_true/num_sequence_keys, ")"
		print "Clusters true: ", sum_clusters_true, " / ", num_sequence_keys, "(", 1.0*sum_clusters_true/num_sequence_keys, ")"
		print "Total true: ", sum_audit_true, " / ", num_sequence_keys, "(", 1.0*sum_audit_true/num_sequence_keys, ")"

		complete_seqs.append(sequence_key) if audit_flag else incomplete_seqs.append(sequence_key)

	return {"complete":complete_seqs, "incomplete":incomplete_seqs}

def audit_sequence(sequence_key):
	mcm_flag, probability_flag = check_mcm(sequence_key)
	convergence_flag = check_convergence(sequence_key)
	clusters_flag = check_clusters(sequence_key)
	
	return (mcm_flag, probability_flag, convergence_flag, clusters_flag)


def check_mcm(sequence_key):

	#kdrew: check to see if 5 entries uploaded into mcm
	new_scop = "1.75"
	new_timestamp = "2010-06-18"
	old_scop = "1.69"
	

	mcmdata_new = session.query(McmData).filter(and_(
		McmData.sequence_key == sequence_key,
		McmData.scop == new_scop,
		McmData.timestamp >= new_timestamp
		)).all()
	new_prob_sum = 0
	for data in mcmdata_new:
		new_prob_sum += data.probability
	avg_new_prob = new_prob_sum/len(mcmdata_new) if len(mcmdata_new) else 0

	mcmdata_old = session.query(McmData).filter(and_(
		McmData.sequence_key == sequence_key,
		McmData.scop == old_scop
		)).all()
	old_prob_sum = 0
	for data in mcmdata_old:
		old_prob_sum += data.probability
	avg_old_prob = old_prob_sum/len(mcmdata_old) if len(mcmdata_old) else 0

	return True if len(mcmdata_new) >= 5 else False, True if avg_new_prob >= avg_old_prob else False

#kdrew: check for rosetta_convergence
def check_convergence(sequence_key):
	convergence_data = session.query(RosettaConvergence,FilesystemOutfile).join(RosettaConvergence.outfile).filter(
		FilesystemOutfile.sequence_key == sequence_key,
		).all()

	#print "convergence_data length: ", len(convergence_data)
	return True if len(convergence_data) == 1 else False

#kdrew: check for rosetta_cluster
def check_clusters(sequence_key):

	cluster_data = session.query(RosettaConvergence,FilesystemOutfile,RosettaCluster).join(RosettaConvergence.outfile,RosettaConvergence.clusters).filter(
		FilesystemOutfile.sequence_key == sequence_key,
		).all()

	#print "cluster_data length: ", len(cluster_data)
	return True if len(cluster_data) >= 25 else False

#kdrew: check for decoys uploaded into structure

if __name__ == "__main__":
	main()
