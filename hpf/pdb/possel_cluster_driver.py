#!/usr/bin/env python

## A driver for the serialization of clustering tasks over families. 
## A "task" is a family, for which clustering analysis is done for
## a set of "sites" (a list of PDB atom numbers)
## 
## Auth: Duncan Penfold-Brown, 9/14/2011

from hpf.utilities import consume
from hpf.processing import processor
from hpf.pdb.cluster import BioPDBStruct
from hpf.hddb.db import Session, ScopedSession, Family, Protein, Domain, Structure, PDBSeqRes

## Globals

usage_str = """ """

samples = 1000
results_file = 'positive_selection_clustering.pkl'


## Option setting and parsing (EXPAND LATER)


## Support functions

def get_possel_sites(family, probability_cutoff=0.0):
# Takes a family and an optional probability cutoff (default 0)
# Returns a list of alignment columns representing pos. selection sites for given family
    if family.selection == None:
        raise Exception("Family {0} has no sites of positive selection".format(family.id))
    sites = []
    for site in family.selection:
        if site.probability >= probability_cutoff:
            sites.append(site.column)
    return sites

def _get_pdb_struct(domain, ):
# Returns a triple (pdb id, chain, structure ORM object) for a PDB (/known struct) domain
    session = ScopedSession()

    # Check for objects needed to query for pdb struct.
    if domain.parent_id == None or domain.sccs == None:
        raise Exception("Domain {0} has no sccs record (cannot find PDB)".format(domain.id))

    # Query to get the PDBSeqRes record matching the domain, return structure info
    psr = session.query(PDBSeqRes).filter_by(sequence_key=domain.parent_id[3:], chain=domain.sccs.chain).first()
    if not psr:
        raise Exception("No PDB record found for domain {0} (seq {1}, chain {2})".format(domain.id, domain.parent_id[3:], domain.sccs.chain))
    elif psr.structure == None:
        raise Exception("Domain {0} PDBSeqRes entry has no structure".format(domain.id))
    return (psr.pdb.pdbId, psr.chain, psr.structure)

def get_domain_struct(domain):
# Returns a tuple of (PDB ID, Chain, and Struct db ORM object) from the given domain
# Throws Exception if no structure available
    if domain.domain_type in ("psiblast", "fold_recognition"):
        (pdb_id, pdb_chain, struct) = _get_pdb_struct(domain)
    elif domain.domain_type in ("msa", "pfam", "unassigned", "user_defined"):
        raise Exception("Currently not considering non-PDB structure-mapped domains")
    else:
        raise Exception("Domain type {1} unknown".format(domain.id, domain.domain_type))
    return (pdb_id, pdb_chain, struct)

def get_domain_sites(family, protein, domain, sites):
# Convert alignment column sites for a family into domain-specific pdb atom num sites.
# NOTE: This currently only is tested for real PDB structures (not decoys)
    from hpf.pdb.mapper import AlignmentToPDBMapper
    pdbmap = AlignmentToPDBMapper(family, protein, domain)
    local_sites = []
    for al_site in sites:
        try:
            local_sites.append(pdbmap.alignment_pdbatom_map[al_site])
        except KeyError:
            #DEBUG 
            print "Alignment column {0} does not exist in PDB record. Ignoring..".format(al_site)
            continue
    return local_sites


## Task functionality

def cluster_driver(family_id):
    session = Session()
    family = session.query(Family).get(family_id)
    if family == None:
        raise Exception("Family {0} could not be fetched from the database".format(family_id))
    
    # Get sites for family (+Sel and TODO: firedb)
    ps_sites = get_possel_sites(family)
    #DEBUG
    print "Start cluster analysis for family {0}".format(family.id)
    print "Family {0} sites of +sel: ".format(family.id), ps_sites
    
    # Get repr protein from family (first one)
    protein = family.proteins[0]

    # Attempt clustering on structures for all domains in protein
    for domain in protein.domains:
        try:
            (pdb_id, pdb_chain, struct) = get_domain_struct(domain)
        except Exception as e:
            print e
            print "Domain {0} (Family {1}, Protein {2}) has no valid structure. Skipping..".format(domain.id, family.id, protein.id)
            continue
        
        # Get domain-specific sites
        domain_sites = get_domain_sites(family, protein, domain, ps_sites)
        #DEBUG
        print "Family {0}, Protein {1}, Domain {2} local sites: ".format(family.id,protein.id,domain.id), domain_sites
        
        # Create BioPDBStruct object to call cluster analysis on
        pdb_struct = BioPDBStruct(pdb_id, pdb_chain, debug=True)
        cluster_id_str = "Family {0}, Protein {1}, Domain {2}, Structure {3}".format(family.id, protein.id, domain.id, struct.id)
        try:
            pdb_struct.cluster_analysis(domain_sites, sample_size=samples, store_file=results_file, tag=cluster_id_str, report=True)
        except Exception as e:
            print e
            print "Can not complete clustering analysis on domain {0}. Skipping..".format(domain.id)
            continue
    
    #DEBUG
    print "Clustering for domain/structs in Family {0} (Protein {1}) complete".format(family.id, protein.id)

def tasks():
    # Task -> a family id to process clustering on
    # Get list of families to cluster (GENERAL version would get list of proteins)
    session = Session()
    #TEST: limit query for testing
    #families = session.query(Family).filter_by(manually_curated=0).limit(5).all()
    families = session.query(Family).filter_by(manually_curated=0).all()
    family_ids = []
    for fam in families:
       family_ids.append(fam.id)
    if family_ids == []:
       raise Exception("No families to cluster were retrieved from the db")
    
    #TEST - manually populate family test set
    #family_ids = [19187, 18883]
    
    return family_ids


def main():
    #TODO: Create a parser for cmdline options
    # Create a processing pool (via hpf.processing) and make tasks to serve to the cluster_driver
    pool = processor()
    pool.make_tasks(tasks, )
    consume(pool.run(cluster_driver))

if __name__ == "__main__":
    main()



