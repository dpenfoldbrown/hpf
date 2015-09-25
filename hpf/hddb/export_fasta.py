#!/usr/bin/env python

## export_fasta
#
# A script to export fasta sequences from the local (HPF) database with local
# header information. Exports into given outfile. Exports all protein sequences 
# from  a list of given experiments
# To format exported fasta for blast use, use formatdb (ncbi)
# dpb 10/24/2011

import sys
from hpf.hddb.db import Session, Experiment, Protein, Sequence

# Export file sequence header format
header = '>hpf|sequence_key|protein_key|experiment_key|experiment_name'

# Globals (experiment list and outfile)
experiments = [1186,]
outfile = 'yeast_1186.hpf.fasta'

def export_fasta(outhandle, experiment):
    """
    outhandle     - a filehandle for writing sequences to
    experiment    - the hpf.experiment DB id to fetch sequences from
    """
    print "Exporting all sequences from Experiment {0}".format(experiment)
    
    session = Session()
    proteins = session.query(Protein).filter_by(experiment_key=experiment)
    num_seqs = proteins.count()

    i = 0
    for protein in proteins:
        outhandle.write(">hpf|{0}|{1}|{2}|{3} ({4})\n".format(protein.sequence_key, protein.id, protein.experiment_key, protein.experiment.name, protein.experiment.short_name))
        outhandle.write("{0}\n".format(protein.sequence.sequence))
        sys.stdout.write("{0}     of {1}\r".format(i, num_seqs))
        i += 1
    
    print "{0} sequences from experiment {1} exported".format(i, experiment)

def main():
    with open(outfile, 'w') as handle:
        #handle.write("{0}\n".format(header))
        for experiment in experiments:
            export_fasta(handle, experiment)
    print "Export complete"


if __name__ == "__main__":
    main()
