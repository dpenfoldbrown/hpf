#!/usr/bin/env python

# A set of support functions for general structure functionality
# Auth." Duncan Penfold-Brown, 10/07/2011

from hpf.hddb.db import Session, ScopedSession, Family, Protein, Domain, Structure, PDBSeqRes, McmData

session = None


def get_domain_structure(domain, probability_cutoff=0.8):
# Takes a hpf.hddb.db Domain object (ORM object), returns a hpf.hddb.db Structure (ORM) object
# Written to allow customization of structure fetching functionality
    # Split the domain get_struct call to PDB or Decoy.
    if domain.domain_type in ("psiblast", "fold_recognition"):
        return _get_pdb_struct(domain)
    elif domain.domain_type in ("msa", "pfam", "unassigned", "user_defined"):
        #return _get_decoy_struct(domain, probability_cutoff)
        print "Currently not considering non-PDB structure-mapped domains. Ignoring domain {0}".format(domain.id)
        return None
    else:
        print "Domain {0} type '{1}' not recognized. Returning None..".format(domain.id, domain.domain_type)
        return None

def _get_pdb_struct(domain, ):
# Returns a DB Structure ORM object
    global session
    if session == None:
        session = Session()

    # Check for objects needed to query for pdb struct.
    if domain.parent_id == None or domain.sccs == None:
        print "Domain {0} has no sccs record (cannot find PDB). Returning None..".format(domain.id)
        return None

    # Query to get the PDBSeqRes record mathing the domain, return structure from PDBSeqRes.
    psr = session.query(PDBSeqRes).filter_by(sequence_key=domain.parent_id[3:], chain=domain.sccs.chain).first()
    if not psr:
        print "No PDB record found for domain {0} (seq {1}, chain {2}). Returning None.".format(domain.id, domain.parent_id[3:], domain.sccs.chain)
        return None
    elif psr.structure == None:
        print "Domain {0} PDBSeqRes entry has no structure. Returning None..".format(domain.id)
    return psr.structure

def _get_decoy_struct(domain, probability_cutoff=0.8):
# Returns a DB Structure ORM object
    global session
    if session == None:
        session = Session()

    # Get decoy with highest MCM value that meets probability cutoff
    mcm = session.query(McmData).filter_by(outfile_key=domain.outfile_key).order_by('probability desc').first()
    if not mcm:
        print "Domain {0} has no decoy structures".format(self.domain.id)
        return None
    elif mcm.probability < probability_cutoff:
        print "Domain {0} structure (prob: {1}) does not meet cutoff {2}".format(domain.id, mcm.probability, probability_cutoff)
        return None
    return mcm.structure

    

