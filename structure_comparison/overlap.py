#
# Overlap functions for domains, astrals, structures, etc.
#
# dpb 2/21/2013
#

def overlap(query_start, query_stop, target_start, target_stop):
    """
    Takes the start and stop points for two regions (converts to ints) and computes the ratio
    of first region (query) that overlaps the second region (target). Returns a float.
    
    Note: the following formula will also compute overlap, but cases are clearer (maybe):
    overlap = [ max(d_stop, a_stop) - min(d_start, a_start) ] - [ abs(d_start - a_start) + abs(d_stop - a_stop) ]
    (negative result => no overlap)
    """

    query_start, query_stop, target_start, target_stop = map(int, (query_start, query_stop, target_start, target_stop))
    query_length = float(query_stop - query_start + 1)

    # No overlap, left or right miss
    if (query_start > target_stop) or (query_stop < target_start):
        overlap = 0.0
    # Query is offset right of target
    elif query_start >= target_start and query_stop >= target_stop:
        overlap = (target_stop - query_start + 1) / query_length
    # Query is offset left of Target
    elif query_start <= target_start and query_stop <= target_stop:
        overlap = (query_stop - target_start + 1) / query_length
    # Query is encapsulated by Target
    elif query_start > target_start and query_stop < target_stop:
        overlap = 1.0
    # Query encapsulates Target
    elif query_start < target_start and query_stop > target_stop:
        overlap = (target_stop - target_start + 1) / query_length
    else:
        raise Exception("Messed up cases in overlap computation")
    if overlap < 0:
        raise ValueError("Overlap of regions {0}-{1} to {2}-{3} is negative".format(query_start, query_stop, target_start, target_stop))
    return float(overlap)

