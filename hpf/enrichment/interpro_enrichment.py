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

## Interpro results file containing results from running background set and enrichment target sequences through interpro
# interpro_results = '/Users/dpb/Documents/superfunc/humanrna/interpro/humanc90s80norna_interpro.out'


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
    sample_counts, enrichment_counts, enrichment_seqs, enrichment_descs = get_counts(enrichment_ids, args.results_file)
    print "Number of Interpro entries in sample: {0}".format(len(sample_counts))
    print "Number of Interpro entries in enrichment list: {0}".format(len(enrichment_counts))

    # Filter enrichment counts and pseudocount enrichment pfams in sample counts
    filter_size = 0
    enrichment_counts = filter_count_by_size(enrichment_counts, filter_size)
    sample_counts = add_psuedo_count(sample_counts, enrichment_counts.keys())

    # Get proportions
    enrichment_proportions = calc_proportion(enrichment_counts)
    sample_proportions = calc_proportion(sample_counts)

    # Calculate enrichment
    interpro_enrichment = calc_enrichment(enrichment_proportions, sample_proportions)

    # Calculate pvalue and bonferroni corrected pvalue
    pvalue_dict = calc_pvalue(enrichment_counts, sample_counts, bonferroni=False)
    corr_pvalue_dict = calc_pvalue(enrichment_counts, sample_counts, bonferroni=True)

    print_sort_pvalue(interpro_enrichment, enrichment_counts, sample_counts, enrichment_seqs, enrichment_descs, pvalue_dict, corr_pvalue_dict)
    print "Complete"


def print_sort_pvalue(enrichment_dict, enrichment_counts, sample_counts, enrichment_seqs, enrichment_descs, pvalue_dict, corr_pvalue_dict):
    # Sort enrichment list by pvalue
    sorted_enrichment = sorted(pvalue_dict.items(), key=lambda k: k[1])
    print "IPR; Enrichment score; # in enrichment set; # in background; IPR desc; p-value; corrected p-value"
    for i in sorted_enrichment:
        print i[0],";", enrichment_dict[i[0]],";", enrichment_counts[i[0]],";", sample_counts[i[0]],";", \
                enrichment_descs[i[0]],";", pvalue_dict[i[0]],";", corr_pvalue_dict[i[0]]

def print_enrichment_list(enrichment_dict, enrichment_counts, sample_counts, enrichment_seqs, enrichment_descs, pvalue_dict, corr_pvalue_dict):
# Print whole list, ordered by enrichment score
    sorted_enrichment = sorted(enrichment_dict.items(), key=lambda k: k[1])
    j = len(sorted_enrichment)
    for i in sorted_enrichment:
        print j,";", i[0],";", i[1],";", enrichment_counts[i[0]],";", sample_counts[i[0]],";", quant_membership(enrichment_seqs[i[0]]), \
                ";", enrichment_descs[i[0]], ";", pvalue_dict[i[0]], ";", corr_pvalue_dict[i[0]]
        j-=1 

def get_counts(enrichment_ids, interpro_outfile):
# Parses through a list of interpro results, building sample and enrichment counts.
# Takes a file of interpro results to count - should contain results from background and enrichment sequences
# Enrichment ids to consider are given
# Returns 3 dicts: sample counts (IPR -> num of occurences),
# enrichment counts (IPR -> num of occurences) for only those query ids in the enrichment list
# enrichment seqs (IPR -> list of seqids that gave that IPR)


    sample_counts = {}
    enrichment_counts = {}
    enrichment_seqs = {}
    enrichment_descriptions = {}
    query_iprs = {}

    # Example results str: 18078925_       EA4EEB033A1433BB        791     HMMPfam PF07679 I-set   166     245     1.2000000000000006E-15  T       25-Oct-2011     IPR013098       Immunoglobulin I-set
    interpro_pattern = r"(?P<query_id>[0-9]+)_?.+[0-9]+\t[0-9]+.+(?P<ipr>IPR[0-9]+|NULL)\t(?P<ipr_desc>[^\t]+)\t*.+$"

    # Build up query -> distinct IPR (set) dict and IPR -> IPR description dict from file
    handle = open(interpro_outfile)
    for line in handle:
        interpro_found = re.match(interpro_pattern, line)
        if not interpro_found:
            raise Exception("Interpro result line '{0}' not matched".format(line))
        query_id = int(interpro_found.group('query_id'))
        ipr = interpro_found.group('ipr')
        try:
            query_iprs[query_id].add(ipr)
        except KeyError:
            query_iprs[query_id] = set()
            query_iprs[query_id].add(ipr)
        if query_id in enrichment_ids and ipr not in enrichment_descriptions.keys():
            enrichment_descriptions[ipr] = interpro_found.group('ipr_desc')
    handle.close()

    # Build IPR counts based on unique set of IPRs per query
    for query_id in query_iprs.keys():
        for ipr in query_iprs[query_id]:
            try:
                sample_counts[ipr] += 1
            except KeyError:
                sample_counts[ipr] = 1
            if query_id in enrichment_ids:
                try:
                    enrichment_counts[ipr] += 1
                    enrichment_seqs[ipr].append(query_id)
                except KeyError:
                    enrichment_counts[ipr] = 1
                    enrichment_seqs[ipr] = [query_id]

    return sample_counts, enrichment_counts, enrichment_seqs, enrichment_descriptions
    
    # ALL IPR RESULTS (not returning only distinct IPRs per query)
    #for line in handle:
    #    interpro_found = re.match(interpro_pattern, line)
    #    if not interpro_found:
    #        raise Exception("Interpro result line '{0}' not matched".format(line))
    #    query_id = int(interpro_found.group('query_id'))
    #    ipr = interpro_found.group('ipr')
    #    try:
    #        sample_counts[ipr] += 1
    #    except KeyError:
    #        sample_counts[ipr] = 1
    #    if query_id in enrichment_ids:
    #        try:
    #            enrichment_counts[ipr] += 1
    #            enrichment_seqs[ipr].append(query_id)
    #        except KeyError:
    #            enrichment_counts[ipr] = 1
    #            enrichment_seqs[ipr] = [query_id]
    #            enrichment_descriptions[ipr] = interpro_found.group('ipr_desc')
    #handle.close()
    #return sample_counts, enrichment_counts, enrichment_seqs, enrichment_descriptions


def ids_from_file(filename):
# Retrieves a list of local (hpf db) sequence IDs from given filename. File must be 
# line-delimited list of sequences only
    handle = open(filename)
    enrichment_ids = []
    for line in handle:
        enrichment_ids.append(int(line.rstrip()))
    handle.close()
    return enrichment_ids

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
