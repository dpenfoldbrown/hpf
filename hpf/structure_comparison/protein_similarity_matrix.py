"""
A script to build protein all-v-all similarity matrices (using the protein-protein
scoring methods found in hpf.structure_comparison.protein_similarity - eg, pairwise
max. average per domain of mammoth zscores).

@auth dpb
@date 12/17/2013
"""

import argparse
from hpf.hddb.db import Session, Protein
from hpf.structure_comparison.protein_similarity import pairwise_max_avg

# Register commandline arguments and grab values
parser = argparse.ArgumentParser(description="Build protein-protein similarity matrix for an HPF experiment")
parser.add_argument("-e", "--experiment", action="store", dest="exp_id", required=True,
    help="HPF Experiment ID (EG: 1186 is Yeast)")
parser.add_argument("-o", "--outfile", action="store", type=argparse.FileType('w'), dest="outhandle", required=True,
    help="File to store generated matrix in (can be large)")
parser.add_argument("-k", "--keyfile", action="store", type=argparse.FileType('w'), dest="keyhandle", required=True,
    help="File to store uniprot+hpf ids in order of their appearance in matrix rows and columns")
args = parser.parse_args()


# Get list of all protein DBOs
session = Session()
proteins = session.query(Protein).filter_by(experiment_key=args.exp_id).all()
if not proteins:
    raise Exception("Error: no proteins for experiment ID {0}".format(args.exp_id))
print "{0} proteins found in experiment {1}".format(len(proteins), args.exp_id)


# Write list info to key file (uniprot ac, hpf seq id, protein id)
args.keyhandle.write("UNIPROT\tHPF_SEQ\tHPF_PROT\n")
for p in proteins:
    if p.ac:
        uniprotac = p.ac.ac
    else:
        uniprotac = " "*5
    args.keyhandle.write("{0}\t{1}\t{2}\n".format(uniprotac, p.sequence_key, p.id))
args.keyhandle.close()


# Iterate symmetrically over list: for each pair i,j, compute i,j similarity,
# store in matrix row for i, write row i to file, del row, move on to next row
for i in range(len(proteins)):
    print "Processing row {0} of {1}".format(i+1, len(proteins))
    row = [0]*len(proteins)
    for j in range(i+1, len(proteins)):
        
        row[j] = pairwise_max_avg(proteins[i], proteins[j])
        
        #try:
        #    row[j] = pairwise_max_avg(proteins[i], proteins[j])
        #except Exception as e:
        #    print "Warning: pairwise max average threw exception: {0}".format(e)
        #    print "\tfor proteins {0}, {1}".format(proteins[i], proteins[j])
        #    row[j] = -1.0
        #    raise e

    for v in row:
        args.outhandle.write("{0}\t".format(v))
    args.outhandle.write("\n")
    args.outhandle.flush()
    print "Completed row {0}".format(i+1)
args.outhandle.close()

print "Writing matrix and keyfile COMPLETE"




