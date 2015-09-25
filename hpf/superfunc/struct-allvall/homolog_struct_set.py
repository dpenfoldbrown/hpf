#!/usr/bin/env python

# A script to read and filter a blast results XML file, and then build a dict of
# { blast query id -> structure object }. The structure object contains blast, 
# protein, and protein structure information, and represents the best blast match
# for the query (where best is considered the best blast match with the best 
# structural coverage, aka 1 astral per match protein domain). 
#
# This defines the structure object, runs the blast parser/filterer, creates the
# dict of structure objects, and stores the dict in a file.
#
# Meant to work with homolog_structure_allvall.py, which does the all v all calcs
# on parts of the query space (hardcoded to do X number of structs v. all)

from hpf.blast.blast_filter import LocalHumanBlastFilter, LocalHumanFiltered
from hpf.hddb.db import Session, ScopedSession, Protein


class HomologComparisonProtein():
# A class to represent protein, structure, and blast (homolog) values for homolog
# structure all-v-all comparison.
    
    def __init__(self, query_id, hit_id, hit_protein, hit_experiment, evalue, bitscore, structures, num_domains):
        self.query_id       = query_id
        self.hit_id         = hit_id
        self.hit_protein    = hit_protein
        self.hit_experiment = hit_experiment
        self.evalue         = evalue
        self.bitscore       = bitscore
        self.structure_list = structures
        self.num_domains    = num_domains

    def __repr__(self,):
        return "Homolog Protein Struct: query {0}, match {1}, eval {2}".format(self.query_id, self.hit_id, self.evalue)


def homolog_struct_records(filtered_records_dict):
# Takes a dict of form query_id -> list[LocalHumanFiltered objs.], where list
# of filtered objects represents blast matches of the query_id query.
# Returns a dictionary of form query_id -> HomologComparisonProtein obj, where
# the HomologComparisonProtein object represents the best blast match with best
# structure coverage from the list of filtered blast matches given.
    count = 0
    no_struct = 0
    homolog_struct_dict = {}
    for query_id in filtered_records_dict:
        count += 1
        
        #DEBUG
        print "Query {0}\t{1}".format(query_id, count)
        
        try:
            query_homolog_rec = get_best_homolog(filtered_records_dict[query_id])
        except LookupError as e:
            print e
            print "Query sequence {0} returned no homologs with good structural coverage".format(query_id)
            no_struct += 1
            continue
        homolog_struct_dict[query_id] = query_homolog_rec
    
    #DEBUG
    print "Total records considered: {0}".format(count)
    print "Records with no structural coverage: {0}".format(no_struct)
    
    return homolog_struct_dict


def get_best_homolog(records):
# Takes a list of Filtered objects. Returns the filtered object with
# best blast score and best structural coverage.
    if not records:
        raise LookupError("Given records list contains no records")
    
    session = ScopedSession()
    
    # Sort records by blast score (record.hit_bitscore), highest first
    records.sort(key=lambda record: record.hit_bitscore, reverse=True)

    best_struct_score = 0.0
    best_record = None
    best_protein = None
    for record in records:
        # Get protein ORM object matching record
        protein = session.query(Protein).get(record.hit_protein_id)
        if not protein:
            print "No protein fetched for query {0} protein ID {1}. Skipping protein..".format(record.query_id, record.hit_protein_id)
            continue
        try:
            protein_score = structure_score(protein)
        except LookupError as l:
            print "Exception {0}. Skipping protein {1}..".format(l, protein.id)
            continue
        
        if protein_score == 1.0:
            best_record = record
            best_protein = protein
            break
        elif protein_score > best_struct_score:
            best_struct_score = protein_score
            best_record = record
            best_protein = protein
    if best_record == None:
        raise LookupError("Could not find a protein with structure")
    
    # Get structure list for best protein
    hit_structures = get_structure_ids(best_protein)
    
    # HomologComparisonProtein prototype:
    # (self, query_id, hit_id, hit_protein, hit_experiment, evalue, bitscore, structures, num_domains)
    return ( HomologComparisonProtein(best_record.query_id, best_record.hit_id, \
                best_record.hit_protein_id, best_record.hit_experiment_id, \
                best_record.hit_evalue, best_record.hit_bitscore, \
                hit_structures, len(best_protein.domains) \
                ) \
            )


def structure_score(protein):
# Takes a protein ORM object. Returns structure coverage score (average of domains'
# best astral overlap. Domain's overlap is 0 if domain has no astral). Returns float.
    if not protein.domains:
        raise LookupError("Protein {0} has no domains".format(protein.id))
    
    total_overlap = 0.0
    for domain in protein.domains:
        if domain.best_astral:
            total_overlap += float(domain.best_astral.overlap)
        # If domain has no astral, do nothing (that domain = 0 overlap)
    return float(total_overlap) / float(len(protein.domains))

def get_structure_ids(protein):
# Takes a protein ORM object. Returns a list of the best astral IDs, one per domain,
# for each domain in the given protein. If a domain has no astral struct, put
# None in results list.
    if not protein.domains:
        raise LookupError("Protein {0} has no domains".format(protein.id))
    
    struct_ids = []
    for domain in protein.domains:
        if domain.best_astral:
            struct_ids.append(domain.best_astral.astral_sid)
        else:
            struct_ids.append(None)
    return struct_ids
        

def store(store_file, dict):
    import cPickle as pickle
    handle = open(store_file, 'w')
    pickle.dump(dict, handle)
    handle.close()


def main():

    # Blast params. TODO: Put everything in config file, if you want to be the coolest.
    #blast_results_file = "/Users/dpb/Documents/superfunc/data/humanrna/human_1176_RNA.c90.s80.blast.xml"
    #blast_results_file = "/Users/dpb/Documents/superfunc/data/humanrna/human_1176_RNA_MORE.c90.s80.blast.xml"
    blast_results_file = "/Users/dpb/Documents/superfunc/data/humanrna/test/new10.blast.xml"
    eval_threshold = 1e-3
    length_threshold   = 0.0
    identity_threshold = 0.0

    #struct_dict_file = 'allhuman_largeblast_query_struct_dict.pkl'
    struct_dict_file = 'test10_id0_len0.pkl'

    # Create blast filter. Parse and filter records
    blast_filter = LocalHumanBlastFilter(eval_cutoff=eval_threshold, \
                        length_cutoff=length_threshold, \
                        identity_cutoff=identity_threshold, \
                        multi_hits=True \
                   )
    blast_records = blast_filter.runBlast(open(blast_results_file))
    filtered_records = blast_filter.filterBlast(blast_records)

    print "Retrieved and filtered blast records. Creating representative structure dict"

    # Get records for structure-structure comparison. Dict, query_id -> HomologComparisonProtein obj.
    #TODO: could change this function to instead just output each struct comparison record to a pickle, one by one
    struct_comparison_records = homolog_struct_records(filtered_records)

    # Store dict of query->protein/structure information to be read and all-v-alled
    store(struct_dict_file, struct_comparison_records)

    print "Creating query->protein structure dict for all-v-all complete"


if __name__ == "__main__":
    main()
