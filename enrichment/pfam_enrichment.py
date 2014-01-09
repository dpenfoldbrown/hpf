#!/usr/bin/env python

# A script to do enrichment analysis of PFAM results.
# Written to parse pfam results files. Builds sample set out of a full set of PFam results.
# Builds enrichment set out of certain target queries in total set (in this case, 700 human RNA protein seqs).
# Copied enrichment functionality from sf_enrichment.py

import re
import sys
import argparse
from hpf.enrichment.enrichment import *
from hpf.hddb.db import Session, ScopedSession, Protein

## Files containing sequence IDs to be used as enrichment sets
# norna_seqfile = '/Users/dpb/bonneau-dev/hpf/trunk/src/projects/berliner_mrna/norna250_seqs.txt'
# novel_humanrna = '/Users/dpb/Documents/superfunc/humanrna_rev2/novel_hpf_seqids.txt'
# humanrnaV2 = '/Users/dpb/Documents/superfunc/humanrna_rev2/humanrnaV2_hpf_seqids.txt'

## Pfam results file containing Pfam results for background set and enrichment set entities (sequences)
# pfam_results = '/Users/dpb/Documents/superfunc/humanrna/pfam/humanc90s80norna_pfam.out'


def main():
    # Create parser and get cmdline arguments
    parser = argparse.ArgumentParser(description="Test argparse module")
    parser.add_argument('-e', '--enrichment_file', metavar="FILE", type=str, action='store', dest='enrichment_file', required=True, 
        help="A line-delim'ed file containing IDs to consider as the enrichment set")
    parser.add_argument('-r', '--results_file', metavar="FILE", type=str, action='store', dest='results_file', required=True, 
        help="The file of results for all background/sample and enrichment IDs")
    #parser.add_argument('-d', '--enrichment_exp', metavar="EXP", type=int, action='store', dest='enrichment_exp', 
    #    help="An HPF DB experiment ID from which to take protein IDs to consider as the enrichment set")
    args = parser.parse_args()
    print "Enrichment set from file: {0}".format(args.enrichment_file)
    print "Results from file: {0}".format(args.results_file)


    # Get enrichment ids to watch for in total set and isolate as enrichment targets
    enrichment_ids = ids_from_file(args.enrichment_file)
    print "Number of enrichment targets: {0}".format(len(enrichment_ids))

    # Get counts
    sample_counts, enrichment_counts, enrichment_seqs, enrichment_descs = get_pfam_counts(enrichment_ids, args.results_file)
    print "Number of PFAM entries in sample (background): {0}".format(len(sample_counts))
    print "Number of PFAM entries in enrichment list: {0}".format(len(enrichment_counts))

    # Filter enrichment counts and pseudocount enrichment pfams in sample counts
    filter_size = 0
    enrichment_counts = filter_count_by_size(enrichment_counts, filter_size)
    sample_counts = add_psuedo_count(sample_counts, enrichment_counts.keys())

    # Get proportions
    enrichment_proportions = calc_proportion(enrichment_counts)
    sample_proportions = calc_proportion(sample_counts)

    # Calculate enrichment
    pfam_enrichment = calc_enrichment(enrichment_proportions, sample_proportions)

    # Calculate pvalue and bonferroni corrected pvalue
    pvalue_dict = calc_pvalue(enrichment_counts, sample_counts, bonferroni=False)
    corr_pvalue_dict = calc_pvalue(enrichment_counts, sample_counts, bonferroni=True)

    print_sort_pvalue(pfam_enrichment, enrichment_counts, sample_counts, enrichment_seqs, enrichment_descs, pvalue_dict, corr_pvalue_dict)
    print "Complete"


def print_sort_pvalue(enrichment_dict, enrichment_counts, sample_counts, enrichment_seqs, enrichment_descs, pvalue_dict, corr_pvalue_dict):
# Sort enrichment list by pvalue
    sorted_enrichment = sorted(pvalue_dict.items(), key=lambda k: k[1])
    print "Pfam; Enrichment score; # in enrichment set; # in background; Pfam name; Pfam desc; p-value; corrected p-value"
    for i in sorted_enrichment:
        print i[0],";", enrichment_dict[i[0]],";", enrichment_counts[i[0]],";", sample_counts[i[0]],";", \
                enrichment_descs[i[0]][0],";", enrichment_descs[i[0]][1],";", pvalue_dict[i[0]],";", corr_pvalue_dict[i[0]]

def print_sort_enrichment(enrichment_dict, enrichment_counts, sample_counts, enrichment_seqs, enrichment_descs, pvalue_dict, corr_pvalue_dict):
# Print whole list, ordered by enrichment score
    sorted_enrichment = sorted(enrichment_dict.items(), key=lambda k: k[1])
    for i in sorted_enrichment:
        print i[0],";", i[1],";", enrichment_counts[i[0]],";", sample_counts[i[0]],";", quant_membership(enrichment_seqs[i[0]]), \
                ";", enrichment_descs[i[0]], ";", pvalue_dict[i[0]], ";", corr_pvalue_dict[i[0]]

def get_pfam_counts(enrichment_ids, pfam_outfile):
# Parses through a list of pfam results, building sample and enrichment counts.
# Takes the filename of a Pfam results file containing Pfam results for background and enrichment seqs
# Enrichment ids to consider are given
# Returns 3 dicts: sample counts (pfam -> num of occurences),
# enrichment counts (pfam -> num of occurences) for only those query ids in the enrichment list
# enrichment seqs (pfam -> list of seqids that gave that pfam)
# enrichment descs (pfam -> description of family)

    sample_counts = {}
    enrichment_counts = {}
    enrichment_seqs = {}
    enrichment_descriptions = {}
    query_pfams = {}

    pfam_pattern = r"(?P<query_id>[0-9]+)_?.+HMMPfam\t(?P<pfam>[a-zA-Z0-9]+)\t(?P<pfam_desc>[^\t]+)\t[0-9]+\t[0-9]+.+(IPR[0-9]+|NULL)\t(?P<ipr_desc>[^\t]+)\t*.+$"

    # Build up query -> distinct Pfam (set) and Pfam -> Pfam Description dict from file
    handle = open(pfam_outfile)
    for line in handle:
        pfam_found = re.match(pfam_pattern, line)
        if not pfam_found:
            raise Exception("PFam result line '{0}' not matched".format(line))
        query_id = int(pfam_found.group('query_id'))
        pfam = pfam_found.group('pfam')
        
        # Add unique Pfams to query set
        try:
            query_pfams[query_id].add(pfam)
        except KeyError:
            query_pfams[query_id] = set()
            query_pfams[query_id].add(pfam)
        
        # If query id is to be reported, store the pfam and IPR descriptions
        if query_id in enrichment_ids and pfam not in enrichment_descriptions.keys():
            enrichment_descriptions[pfam] = (pfam_found.group('pfam_desc'), pfam_found.group('ipr_desc'))
    handle.close()
        
    # Build Pfam counts based on the number of unique Pfams per query    
    for query_id in query_pfams.keys():
        for pfam in query_pfams[query_id]:
            try:
                sample_counts[pfam] += 1
            except KeyError:
                sample_counts[pfam] = 1
            if query_id in enrichment_ids:
                try:
                    enrichment_counts[pfam] += 1
                    enrichment_seqs[pfam].append(query_id)
                except KeyError:
                    enrichment_counts[pfam] = 1
                    enrichment_seqs[pfam] = [query_id]

    return sample_counts, enrichment_counts, enrichment_seqs, enrichment_descriptions

def ids_from_file(filename):
# Retrieves a list of local (hpf db) sequence IDs from given filename. File must be 
# line-delimited list of sequences only
    handle = open(filename)
    ids = []
    for line in handle:
        ids.append(int(line.rstrip()))
    handle.close()
    return ids

def ids_from_db():
# Retrieves a list of local (hpf db) sequence IDs for the sequences you want to calculate enrichment on
# These must be the same ID as in the ID field of pfam/interpro results

    # Get sequence IDs for all humanrna proteins (773 list)
    session = Session()
    proteins = session.query(Protein).filter_by(experiment_key=1177).all()
    protein_seqids = []
    for protein in proteins:
        protein_seqids.append(protein.sequence_key)
    return protein_seqids


if __name__ == "__main__":
    main()
