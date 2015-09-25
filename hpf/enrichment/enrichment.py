#!/usr/bin/env python

# A support file containing all enrichment functions.
# The enrichment process is usually as follows:
#   1) Create an enrichment count, as in: enrichment id -> # of times enrichment id occurs in results
#   2) Create a background (sample count): id -> # of times id occurs in results
#       - Note that your enrichment set should be some set of interest, generally a subset of the background
#   3) Remove low-liers from your enrichment and sample count sets if desired (remove ids with few or no results)
#       - filter_count_by_size(unfiltered_count_dict, threshold)
#   4) Add a pseudocount for the IDs in the enrichment set to the sample set
#       - add_pseudo_count(sample_dict, enrichment_ids, count=1)
#   5) Create proportion dicts of both enrichment and sample sets. 
#       - calc_propotion(count_dict)
#   6) Calculate enrichment scores for enrichment set against sample set.
#       - calc_enrichment(enrichment_proportions, sample_proportions, log=True)
#   7) Calculate pvalues for the enrichment set against the background set
#       - calc_pvalue(enrichment_count_dict, sample_count_dict, bonferroni=False)


from numpy import average
import rpy2.robjects as robj

TINY_NUM = 0.0000000001


def filter_count_by_size(unfiltered_count_dict, threshold):
# Takes a dictionary of counts, returning a new dictionary of only those
# entries with a count >= given threshold
# unfiltered_count_dict of form key -> (int) count
    new_dict = dict()
    for key in unfiltered_count_dict.keys():
        if unfiltered_count_dict[key] >= threshold:
            new_dict[key] = unfiltered_count_dict[key]
    return new_dict


def add_psuedo_count(sample_dict, enrichment_ids, count=1):
# Takes a sample space (background) count dictionary, and adds a given pseudocount
# (default 1) to the sample space entries that are in the given enrichment list. If an
# enrichment entry is not in the sample space, adds enrichment entry to the sample space
# with value 'count'
    for key in enrichment_ids:
        try:
            sample_dict[key] += count
        except KeyError:
            sample_dict[key] = count
    return sample_dict


def calc_proportion(count_dict):
# Takes a dictionary of form key -> # of instances (count)
# Adds the # of instances for all keys, and creates a propotion dictionary of the form 
# key -> % of total count. Returns proportion dictionary.
    total = 0
    proportion_dict = dict()
    for sf in count_dict.keys():
        total += count_dict[sf]
    for sf in count_dict.keys():
        proportion_dict[sf] = float(count_dict[sf])/total
    return proportion_dict


def calc_enrichment(enrichment_proportions, sample_proportions, log=True):
# Calculates the enrichment score for each given enrichment key -> proportion pair.
# Takes a dictionary of enrichment porportions (enrichment id -> proportion) and a
# dictionary of sample space (background) proportions (id -> proportion). Returns
# a dictionary of enrichment scores, enrichment id -> enrichment score
    enrichment_dict = dict()
    for key in enrichment_proportions:
        try:
            enrichment_dict[key] = enrichment_proportions[key] / sample_proportions[key]
        except KeyError:
            #kdrew: this just means the sampling dict does not have the superfamily and therefore it is high enrichment
            #dpb: should not happen as long as pseudocounts for enrichment ids are added to sample counts
            enrichment_dict[key] = enrichment_proportions[key] / TINY_NUM
        if log:
            import math
            enrichment_dict[key] = math.log(enrichment_dict[key])
    return enrichment_dict


def calc_pvalue(counts_dict1, counts_dict2, bonferroni=False):
# Calculates pvalues for each key in the given counts_dict1. Can do bonferroni correction.
# Takes two dictionaries of the form id -> (int) count
# For enrichment, counts_dict1 should be the enrichment counts (# of instances of enrichment ids in results),
# while counts_dict2 should be the sample (or background) counts
# Returns a dict of pvalues for each key in counts_dict1
    pvalue_dict = dict()
    #kdrew: find the total number of counts in the set (dict1) and total set+nonset (dict2)
    set_sum = 0
    total_sum = 0
    bonferroni_count = 0
    for i in counts_dict1:
        set_sum += counts_dict1[i]
        bonferroni_count +=1
    for i in counts_dict2:
        total_sum += counts_dict2[i]

    #print set_sum
    #print total_sum
    #print bonferroni_count
    
    for key in counts_dict1:
        #kdrew: set up r matrix, assumes set (dict1) is a subset of total (dict2)
        sf_set = counts_dict1[key]
        sf_nonset = counts_dict2[key] - counts_dict1[key]
        nonsf_set = set_sum - counts_dict1[key]
        nonsf_nonset = total_sum - sf_nonset
        v = robj.IntVector([sf_set, nonsf_set, sf_nonset, nonsf_nonset])
        m = robj.r['matrix'](v,nrow=2)
        ftest = robj.r['fisher.test']
        ft = ftest(m, alternative="greater")
        if bonferroni:
            pvalue_dict[key] = ft[0][0]*bonferroni_count
        else:
            pvalue_dict[key] = ft[0][0]
    return pvalue_dict


