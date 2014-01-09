##
## Functions for representing Domains with actual structures.
#
## For known structures, returns structure keys of astral structures that map onto
## given domain by way of domain -> pdb -> astral overlap mapping. For known structures,
## the returned structure keys are meant to collectively represent the domain.
#
## For "unknown" structures, returns a list of the structure keys of the denovo structures
## from the best 5 MCM score denovo structures. In this case, the structures are meant to
## individually represent the domain.
#
## dpb 5/01/2013

import warnings

# Cutoffs for astral finding
ASTRAL_MCM_CUTOFF = 0.8
ASTRAL2DOMAIN_OVERLAP_CUTOFF = 0.8
ASTRAL2ASTRAL_OVERLAP_CUTOFF = 0.1


def structure_representation(domain):
    """
    Parameters:
        domain  -   hpf.hddb.db.Domain object
    Return:
        [int]   -   list of hpf.hddb.Structure keys (IDs)
    """
    if domain.domain_type in ("psiblast", "fold_recognition"):

        if domain.astral_domain_overlap:
            # Remove astrals with overlap less than cutoff
            #thresholded_astral_overlaps = [a for a in domain.astral_domain_overlap if a.overlap >= ASTRAL2DOMAIN_OVERLAP_CUTOFF]
            thresholded_astral_overlaps  = list()
            for a in domain.astral_domain_overlap:
                if float(a.overlap) >= ASTRAL2DOMAIN_OVERLAP_CUTOFF:
                    thresholded_astral_overlaps.append(a)

                    #print "DEBUG: Adding {0} to thresholded astrals".format(a)
                    #print "DEBUG: overlap: {0}, >= {1}: {2}".format(a.overlap, ASTRAL2DOMAIN_OVERLAP_CUTOFF, a.overlap >= ASTRAL2DOMAIN_OVERLAP_CUTOFF)
    
            # Remove Astrals with same structure key (just in case)
            astral_struct_dict = dict()
            for a in thresholded_astral_overlaps:
                astral_struct_dict[a.astral.structure_key] = a
            clean_astral_overlaps = astral_struct_dict.values()
    
            if (len(clean_astral_overlaps) < 2):
                return [a.astral.structure_key for a in clean_astral_overlaps]
            
            #print "\tDEBUG: Clean astrals: ", clean_astral_overlaps
            
            # Remove astrals with same astral_start/astral_stop pair as those already in the set (further cleaning step)
            start_stop_pairs = []
            nonredun_astral_overlaps = []
            for a in clean_astral_overlaps:
                if (a.astral_start, a.astral_stop) not in start_stop_pairs:
                    start_stop_pairs.append((a.astral_start, a.astral_stop))
                    nonredun_astral_overlaps.append(a)
            clean_astral_overlaps = nonredun_astral_overlaps

            # Build up graph of thresholded astrals. 
            ## Node ID is astral structure key 
            ## Score is # of domain residues covered by astral (astral.overlap * astral.length)
            import networkx
            from itertools import combinations
            from hpf.graph import max_score_path, path_score
            from hpf.structure_comparison.overlap import overlap
            dag = networkx.DiGraph()
            for a in clean_astral_overlaps:
                dag.add_node(a.astral.structure_key, score=(a.overlap * (a.astral_stop - a.astral_start + 1)))
    
            #print "\tDEBUG: Graph nodes: ", dag.nodes()
    
            ## Build up graph edges. Nodes are connected if the are "complimentary" (non- or nearly non-overlapping in domain space)
            for (a, b) in combinations(clean_astral_overlaps, 2):
                a2a_overlap = overlap(a.astral_start, a.astral_stop, b.astral_start, b.astral_stop)
                if a2a_overlap < ASTRAL2ASTRAL_OVERLAP_CUTOFF:
                    dag.add_edge(a.astral.structure_key, b.astral.structure_key)
    
            #print "\tDEBUG: Graph edges: ", dag.edges()
   
            # Check for overlarge graph before starting recursive methods
            if len(dag.edges()) > 50:
                warnings.warn("Too many redundant nodes in graph for domain {0}, skipping..".format(domain))
                return []

            # Find max scoring path from all nodes, keeping max (could be written better to store as graph was built. Oh well)
            ## NOTE: a path is a list of nodes, and a node is just an astral ID (when fetched via graph.get_nodes()
            max_score = -1
            max_path = None
            for node in dag.nodes():
                path = max_score_path(dag, node)
                score = path_score(dag, path)
    
                #print "DEBUG: Path {0}, Score {1}".format(path, score)
    
                if score > max_score:
                    max_score = score
                    max_path = path
            return max_path if max_path else []

    elif domain.domain_type in ("pfam", "msa", "unassigned"):
        if domain.mcmdata:
            # Remove MCM entries with the same struct key, keeping that with the better probability
            mcm_dict = dict()
            for m in domain.mcmdata:
                if m.structure_key in mcm_dict.keys():
                    if m.probability > mcm_dict[m.structure_key].probability:
                        mcm_dict[m.structure_key] = m
                else:
                    mcm_dict[m.structure_key] = m
    
            # Sort the unique MCM entries by probability (high to low), keeping the highest 5
            mcms = mcm_dict.values()
            if (not mcms):
                return []
            mcms = sorted(mcm_dict.values(), key=lambda k: float(k.probability), reverse=True)[:5]
    
            # Return structure keys for highest 5 non-duplicate MCMs
            return [m.structure_key for m in mcms]
    return []

