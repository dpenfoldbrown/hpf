#
# Init library for structure comparison module
#
# dpb 2/21/2013
#

def get_astrals(domain, session):
    """
    Takes a domain ORM object and a DB session.
    Returns a list of all unique astrals (DBOs) matching ALL given domain's pdbSeqRes records
    """
    from hpf.hddb.db import Astral
    all_astrals = []
    for psr in domain.pdbseqs:
        psr_astrals = session.query(Astral).filter_by(pdbid=psr.pdb.pdbId).filter(Astral.chain.like(psr.chain + '%')).all()
        for ast in psr_astrals:
            if ast not in all_astrals:
                all_astrals.append(ast)
    return all_astrals

def get_astral_startstop(astral):
    """
    Takes an Astral DBO (from hpf orm hpf.hddb.db), attempts to parse start stop from chain value.
    If chain has no range (indicating astral covers whole PDB), get range from astral sequence length
    Returns a tuple of ints, (start, stop)
    """
    import re
    range_chain_pat = r"[0-9A-Z]{1}:(?P<negative>-?)(?P<start>[0-9]+)[A-Z]?-(?P<stop>[0-9]+)[A-Z]?.*"
    empty_chain_pat = r"[0-9A-Z]{1}:$"
    chain_found = re.match(range_chain_pat, astral.chain)
    if chain_found:
        if chain_found.group('negative'):
            return(0, int(chain_found.group('stop')))
        else:
            return (int(chain_found.group('start')), int(chain_found.group('stop')))
    elif re.match(empty_chain_pat, astral.chain):
        return (1, len(astral.sequence.sequence))
    raise Exception("Astral chain string {0} does not match known formats".format(astral.chain))

def parse_astral_chain(chain):
    """
    Parse the astral chain string (eg.s 'B:', 'A:1-23', 'C:12-21,30-120'), returns just chain letter
    """
    import re
    chain = chain.split(':', 1)[0]
    if not re.match(r"[0-9A-Z]{1}", chain):
        raise Exception("Parsed astral chain '{0}' is not a single letter or number".format(chain))
    return chain
