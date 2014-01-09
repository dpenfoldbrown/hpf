#!/usr/bin/env python

# A script to parse XML blast results, filter those results with hpf.blast.blast_filter,
# then retrieve the best blast match per query sequence that has a confident domain_sccs
# entry. For enrichment, output a set (list) of protein sequence keys representing best 
# blast matches.

from hpf.hddb.db import ScopedSession, Protein
from hpf.blast.blast_filter import HumanRNABlastFilter, LocalHumanFiltered


def get_best_record_key(records):
# Takes records, a list of LocalHumanFiltered objects.
# Returns the protein sequence key of the "best" filtered record, where best
# is (th.) the best blast hit with the highest confident sccs entry
# ACTUAL return is a triple: for the best protein, (query id, seq key, prot key, experiment key)

    if not records:
        raise Exception("Given records list contains no records")

    #DEBUG
    for rec in records:
        print "\tHit: {0}, {1}, {2}, eval: {3} bitscore {4}".format(rec.hit_id, rec.hit_protein_id, rec.hit_experiment_id, rec.hit_evalue, rec.hit_bitscore)

    session = ScopedSession()
    
    # Sort records by blast score (record.hit_bitscore), highest first
    records.sort(key=lambda record: record.hit_bitscore, reverse=True)

    # record is LocalHumanFiltered obj: query_id, hit_id, hit_protein_id, hit_experiment_id, hit_evalue, hit_bitscore
    best_sccs_conf = 0.0
    best_record = None
    for record in records:
        
        # Get protein ORM object matching record
        protein = session.query(Protein).get(record.hit_protein_id)
        if not protein:
            raise Exception("No protein fetched for record protein ID {0}".format(record.hit_protein_id))
        try:
            prot_best_sccs = _get_best_sccs(protein)
        except Exception as e:
            print e
            continue
        
        if float(prot_best_sccs.confidence) == 1.0:
            return (record.query_id, protein.sequence_key, protein.id, protein.experiment_key)
        elif float(prot_best_sccs.confidence) > best_sccs_conf:
            best_sccs_conf = float(prot_best_sccs.confidence)
            best_record = record
    if best_record == None:
        raise Exception("Could not find a best protein sequence")
    
    hit_id, hit_protein_id, hit_experiment_id = map(int, (best_record.hit_id, best_record.hit_protein_id, best_record.hit_experiment_id))
    return (record.query_id, hit_id, hit_protein_id, hit_experiment_id)


def _get_best_sccs(protein):
# Takes a Protein ORM object (hpf.hddb.db), returns the sccs object attached to given
# protein with best confidence value
    if not protein.sccs:
        raise Exception("Protein {0} has no sccs records".format(protein.id))
    
    best_conf = 0.0
    best_sccs = None
    for sccs in protein.sccs:
        if float(sccs.confidence) == 1.0:
            return sccs
        elif float(sccs.confidence) >= best_conf:
            best_conf = float(sccs.confidence)
            best_sccs = sccs
    if not best_sccs:
        raise Exception("Protein {0}'s sccs entries not valid (all < 0)".format(protein.id))
    return best_sccs

def list_to_file(list, file):
    outhandle = open(file, 'w')
    for el in list:
        outhandle.write("{0}\n".format(el))
    outhandle.close()

def list_to_pkl(list, file):
    import cPickle as pickle
    outhandle = open(file, 'w')
    for el in list:
        pickle.dump(el, outhandle)
    outhandle.close()


def main():

    blast_results_file = "/Users/dpb/Documents/superfunc/data/humanrna/human_rna.blast.out.xml"
    eval_threshold = 1e-3
    identity_threshold = 0.5
    length_threshold   = 0.8
    out_file = "human_rna_sccs.blast.enrichment.pkl"

    # Create a blast filter (custom for human RNA dataset)
    blast_filter = HumanRNABlastFilter(eval_cutoff=eval_threshold, \
                        length_cutoff=length_threshold, \
                        identity_cutoff=identity_threshold, \
                        multi_hits=True \
                   )

    # Call blast filter for records, then filter records (filter records return dict)
    blast_records = blast_filter.runBlast(open(blast_results_file))
    filtered_records = blast_filter.filterBlast(blast_records)

    #DEBUG
    #print "Filtered records:"
    #blast_filter.print_filtered_records(filtered_records)

    # For each query sequence, get the protein sequence key of the best (highest blast score & best sccs entry) hit
    enrichment_entries = []
    for query_id in filtered_records:
        #DEBUG
        print "Processing blast hits for query sequence {0}".format(query_id)
        try:
            query_id, prot_seq, prot_id, experiment = get_best_record_key(filtered_records[query_id])
        except Exception as e:
            print e
            print "Query sequence {0} returned no protein sequence key for enrichment".format(query_id)
            continue
        enrichment_entries.append( (query_id, prot_seq, prot_id, experiment) )
    
    if enrichment_entries:
        list_to_pkl(enrichment_entries, out_file)
        #list_to_file(enrichment_seqids, out_file)
    else:
        print "No homologous sequences for enrichment found in blast results '{0}'".format(blast_results_file)
    print "Completed processing enrichment homologs for blast results in {0}".format(blast_results_file)

if __name__ == "__main__":
    main()


