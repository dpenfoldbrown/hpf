"""
Protein vs protein similarity functions
"""

from hpf.structure_comparison.parse_structcmp_db import parse_all, StructPair
from hpf.structure_comparison.structure_mammoth import structure_mammoth, run_mammoth

# GLOBAL
global_struct_sets = {}

# FUNCTIONS
def get_structset(protein, use_global=False):
    """
    Gets the representative structure set of a given protein (all structures in nested format
    representing the protein). Returns a nested list of structure keys.
    Nesting is as follows: For protein domains of known type, a number (1+) struct keys will be added
    flatly to the protein's structure set. For protein domains of unknown (decoy) type, a "nested" list 
    of structure keys will be added to the protein's structure set.
    Parameters:
        protein -   protein DBO (hpd.hddb.db.Protein)
    Return:
        list    -   List of elements, where element is structure key or [structure keys]
    """
    global global_struct_sets
    if use_global:
        if protein.id in global_struct_sets:
            #print "Using cached protein structure set for protein {0}".format(protein.id)
            return global_struct_sets[protein.id]
        #else:
        #    print "No cached protein structure set for protein {0}".format(protein.id)
        
    from hpf.structure_comparison.domain_structure_representation import *
    struct_set = []
    for domain in protein.domains:
        domain_structure_set = structure_representation(domain)
        if not domain_structure_set:
            # Don't add to protein struct set if empty
            continue
        elif domain.known_type:
            struct_set += domain_structure_set
        else:
            struct_set.append(domain_structure_set)
    
    if use_global:
        global_struct_sets[protein.id] = struct_set

    return struct_set


def is_known(structval):
    """
    Returns true if structval from structset is a single struct key value, from a known domain
    (so, false if the structval is a list of struct keys, as from an unknown domain)
    """
    return not isinstance(structval, list)


def mem_structure_similarity(source, target, source_type, target_type):
    """In-memory version of structure similarity"""
    sp = StructPair(source=str(source), target=str(target))
    if sp in global_sm_db:
        #print "Found source {0}, target {1} in memory".format(source, target)
        return float(global_sm_db[sp])

    sp = StructPair(source=str(target), target=str(source))
    if sp in global_sm_db:
        #print "Found source {0}, target {1} in memory".format(target, source)
        return float(global_sm_db[sp])

    # If not in mem, run mammoth and return zscore
    #print "Source {0}, target {1} not in memory dict, mammothing".format(source, target)
    mammoth_result = run_mammoth(source, target, source_type, target_type)
    if not mammoth_result:
        raise Exception("RunMammoth failed to return score for structures {0}, {0}".format(source, target))
    
    return float(mammoth_result.zscore)


def structure_similarity(source, target, source_type, target_type):
    """
    Returns the Mammoth similarity score of two structures.
    Parameters:
        source  -   (int) structure key of source structure
        target  -   (int) structure key of target structure
        source_type - (str) the "type" of source structure (see hpf.structure DB table for types) 
        target_type - (str) the type of target structure (%)
    Return:
        float   -   the Mammoth "score" (zscore) of the two structures
    """
    
    # DON'T DO THIS. ENDS UP BEING VERY SLOW IF CALLED MANY TIMES
    # Get structure objects for both struct keys
    #source_dbo = session.query(Structure).get(source)
    #target_dbo = session.query(Structure).get(target)
    #if not (source and target):
    #    raise Exception("One of source {0} or target {1} structures not found in DB".format(source, target))

    # Get mammoth record (checks DB, retrieves if exists. Else runs mammoth and returns)
    sm = structure_mammoth(source,
                      target,
                      source_type,
                      target_type,
                      dbstore=True,
                      cleanup=True,
                      debug=False,
                      table_destination="part01")
    if not sm:
        raise Exception("Structure Mammoth failed to return a mammoth score for structures {0}, {1}".format(source, target))

    return float(sm.zscore)

def pairwise_max_avg(psource, ptarget):
    """
    Basic pairwise max average. Calculates the similarity between two proteins.
    Parameters:
        psource -   protein DBO (hpf.hddb.db.Protein) of "source" comparison protein
        ptarget -   protein DBO (hpf.hddb.db.Protein) of "target" comparison protein
    Return:
        float   -   "score" of calculated protein v protein similarity     
    """
    #DEBUG
    #print "Source: {0}\t Target: {1}".format(psource, ptarget)

    forward_scores = unidirectional_pairwise_max_avg(psource, ptarget)
    backward_scores = unidirectional_pairwise_max_avg(ptarget, psource)
    scores = forward_scores + backward_scores 
    
    if len(scores) < 1:
        #print "Warning: No scores returned to pairwise_max_average"
        return 0.0

    # Sum and average
    score = sum(scores) / len(scores)

    #DEBUG
    #print "Forward scores:  {0}".format(forward_scores)
    #print "Backward scores: {0}".format(backward_scores)

    return score

def unidirectional_pairwise_max_avg(psource, ptarget):
    """
    Basic pairwise max average. Calculates the similarity between two proteins.
    Parameters:
        psource -   protein DBO (hpf.hddb.db.Protein) of "source" comparison protein
        ptarget -   protein DBO (hpf.hddb.db.Protein) of "target" comparison protein
    Return:
        list[float] -   list of the maximum pairwise scores from source to target protein   
    """
    
    source_structset = get_structset(psource, use_global=True)
    target_structset = get_structset(ptarget, use_global=True)
    
    #DEBUG
    #print "Source protein struct set: {0}".format(source_structset)
    #print "Target protein struct set: {0}".format(target_structset)

    max_pair_scores = []
    for sstruct in source_structset:
        source_scores = []
        for tstruct in target_structset:
            source_score = 0
            if is_known(sstruct):   
                if is_known(tstruct):       # BOTH KNOWN
                    source_score = SIMILARITY_METHOD(sstruct, tstruct, "astral", "astral")
                else:                       # Source KNOWN, Target UNKNOWN
                    max_pair_score = 0
                    for tstruct_decoy in tstruct:
                        score = SIMILARITY_METHOD(sstruct, tstruct_decoy, "astral", "decoy")
                        if score > max_pair_score:
                            max_pair_score = score
                    source_score = max_pair_score
            else:
                if is_known(tstruct):       # Source UNKNOWN, Target KNOWN
                    max_pair_score = 0
                    for sstruct_decoy in sstruct:
                        score = SIMILARITY_METHOD(sstruct_decoy, tstruct, "decoy", "astral")
                        if score > max_pair_score:
                            max_pair_score = score
                    source_score = max_pair_score
                else:                       # BOTH UNKNOWN
                    max_pair_score = 0
                    for sstruct_decoy in sstruct:
                        for tstruct_decoy in tstruct:
                            score = SIMILARITY_METHOD(sstruct_decoy, tstruct_decoy, "decoy", "decoy")
                            if score > max_pair_score:
                                max_pair_score = score
                    source_score = max_pair_score
            source_scores.append(source_score)

        # Take the max of all the scores for this source structure
        if len(source_scores) > 0:
            max_source_score = max(source_scores)
            max_pair_scores.append(max_source_score)

        #DEBUG
        #print "Source {0}\tMax Score {1}".format(sstruct, max_source_score)        
    
    return max_pair_scores


#SIMILARITY_METHOD = structure_similarity
SIMILARITY_METHOD = mem_structure_similarity

if SIMILARITY_METHOD == mem_structure_similarity:
    print "IMPORTANT! Parsing structure mammoth DBs into memory. This will take some time..."
    global_sm_db = parse_all()
    print "Completed parsing structure mammoth DBs"

if __name__ == "__main__":
    print "Main not implemented"
    pass
