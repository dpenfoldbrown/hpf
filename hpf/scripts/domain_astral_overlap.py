#!/usr/bin/env python

# A script to score domains with known structure (PDBs) in regard to how much they overlap
# with astral structures that make up parts of the domain's known structure.
# IE: for each domain, get known struct pdb. for the pdb, compare region of pdb domain maps 
#   to region of pdb different astral structures map to. Score based on percent of domain
#   that overlaps with each astral.

# Can store overlap to pickle, and can add to database hpf.domain_astral_overlap
# Based on experiment ID

import re
import cPickle as pickle
from sqlalchemy.exc import IntegrityError
from hpf.hddb.db import Session, ScopedSession, Protein, Domain, Astral, PDBSeqRes, DomainAstralOverlap


def domain_to_astral(experiment_id, store_file=None, db_import=False):

    # Open store pickle if given (to pass to store function later)
    if store_file:
        store_handle = open(store_file, 'a')
    
    # Create session and retrieve domains with known structure for given experiment
    session = Session()
    domains = session.query(Domain).join(Protein).filter_by(experiment_key=experiment_id).filter(Domain.domain_type.in_(['psiblast', 'fold_recognition'])).all()
    print "Considering {0} domains to match to astrals".format(len(domains))

    # For each domain, get PDB info and domain start/stop points, find list of matching Astrals, and do overlap calculation
    count = 0
    no_astral = 0
    num_domains = len(domains)
    
    for domain in domains:
        count += 1
        print "{0}  of  {1}".format(count, num_domains)
        
        # Set domain start/stop values
        domain_start = domain.region.parent_start
        domain_stop  = domain.region.parent_stop

        # Get all astral entries for all domain's pdb IDs
        astrals = get_astrals(domain, session)
        if not astrals:
            print "No astrals found for domain {0}".format(domain.id)
            no_astral += 1
            continue

        # For each astral, get PDB start/stop points, calc domain overlap, print and store
        for astral in astrals:
            try:
                (astral_start, astral_stop) = parse_astral_startstop(astral.chain)
                overlap_ratio = overlap(domain_start, domain_stop, astral_start, astral_stop)
            except ValueError:
                print "Negative overlap for domain {0} ({1}-{2}), Astral {3} (PDB {4}{5})".format(domain.id, domain_start, domain_stop, astral.sid, astral.pdbid, astral.chain)
                raise
            except:
                print "Error calculating overlap for  domain {0} ({1}-{2}), Astral {3} (PBD {4}{5})".format(domain.id, domain_start, domain_stop, astral.sid, astral.pdbid, astral.chain)
                raise
            
            # Print overlap details and store in pickle file
            print "Domain {0} ({1}-{2}), Astral {3} (PDB {4}{5}), Overlap {6}".format(domain.id, domain_start, domain_stop, astral.sid, astral.pdbid, astral.chain, overlap_ratio)
            if store_file:
                store_domainastral_overlap(domain.id, domain_start, domain_stop, astral.id, astral.sid, astral_start, astral_stop, astral.pdbid, astral.chain, overlap_ratio, store_handle)
            
            # Push record into database (unless domain and astral do not overlap)
            if db_import and overlap_ratio > 0.0:
                chain = parse_astral_chain(astral.chain)
                try:
                    db_add(domain.id, astral.id, astral.sid, domain_start, domain_stop, astral_start, astral_stop, astral.pdbid, chain, overlap_ratio, session=session)
                except Exception as e:
                    print e
                    print "Adding domain {0}, astral {1} entry to database failed".format(domain.id, astral.id)
                    raise

    # Finish and close
    if store_file: store_handle.close()
    print "Domain-Astral overlap calculation complete for experiment {0}".format(experiment_id)
    print "{0} domains returned no astral structures".format(no_astral)
    if store_file: print "Results stored in {0}".format(store_file)


def get_astrals(domain, session):
# Takes a domain ORM object and a DB session (see note on get_pdbid function).
# Returns a list of all unique astrals matching ALL given domain's pdbSeqRes records
    all_astrals = []
    for psr in domain.pdbseqs:
        psr_astrals = session.query(Astral).filter_by(pdbid=psr.pdb.pdbId).filter(Astral.chain.like(psr.chain + '%')).all()
        for ast in psr_astrals:
            if ast not in all_astrals:
                all_astrals.append(ast)
    return all_astrals

def parse_astral_startstop(chain):
# Takes an Astral object's chain, parse and return tuple (start, stop).
# Returns (None, None) if chain has no start stop values (indicating astral matches the whole PDB)
    range_chain_pat = r"[0-9A-Z]{1}:(?P<negative>-?)(?P<start>[0-9]+)[A-Z]?-(?P<stop>[0-9]+)[A-Z]?.*"
    empty_chain_pat = r"[0-9A-Z]{1}:$"
    chain_found = re.match(range_chain_pat, chain)
    if chain_found:
        if chain_found.group('negative'):
            return('0', chain_found.group('stop'))
        else:
            return (chain_found.group('start'), chain_found.group('stop'))
    elif re.match(empty_chain_pat, chain):
        return (None, None)
    raise Exception("Astral chain string {0} does not match known formats".format(chain))

def parse_astral_chain(chain):
# Parse the astral chain string (eg.s 'B:', 'A:1-23', 'C:12-21,30-120') to get just the chain letter
    chain = chain.split(':', 1)[0]
    if not re.match(r"[0-9A-Z]{1}", chain):
        raise Exception("Parsed astral chain '{0}' is not a single letter or number".format(chain))
    return chain

def overlap(query_start, query_stop, target_start, target_stop):
# Takes the start and stop points for two regions (converts to ints) and computes the ratio
# of first region that overlaps the second region. Returns a float.
# Note: this function will also compute overlap, but cases are clearer (maybe):
# overlap = [ max(d_stop, a_stop) - min(d_start, a_start) ] - [ abs(d_start - a_start) + abs(d_stop - a_stop) ]
# (negative result => no overlap)

    # Astral ranges are none (Astral overlaps whole PDB): 1.0
    if target_start == None:
        return 1.0
    query_start, query_stop, target_start, target_stop = map(int, (query_start, query_stop, target_start, target_stop))
    dom_length = float(query_stop - query_start + 1)

    # No overlap, left or right miss
    if (query_start > target_stop) or (query_stop < target_start):
        overlap = 0.0
    # Domain is offset right of astral
    elif query_start >= target_start and query_stop >= target_stop: 
        overlap = (target_stop - query_start + 1) / dom_length
    # Domain is offset left of astral
    elif query_start <= target_start and query_stop <= target_stop:
        overlap = (query_stop - target_start + 1) / dom_length
    # Domain is encapsulated by astral
    elif query_start > target_start and query_stop < target_stop:
        overlap = 1.0
    # Domain encapsulates astral
    elif query_start < target_start and query_stop > target_stop:
        overlap = (target_stop - target_start + 1) / dom_length
    else:
        raise Exception("Messed up cases in overlap computation")
    if overlap < 0:
        raise ValueError("Overlap of regions {0}-{1} to {2}-{3} is negative".format(query_start, query_stop, target_start, target_stop))
    return overlap

def store_domainastral_overlap(dom_id, dom_start, dom_stop, astral_id, astral_sid, astral_start, astral_stop, pdb, astral_chain, overlap, store_file_handle):
# Stores all fields in given pickle file
    chain = parse_astral_chain(astral_chain)
    pickle.dump( (dom_id, dom_start, dom_stop, astral_id, astral_sid, astral_start, astral_stop, pdb, chain, overlap), store_file_handle)

def db_add(domain_id, astral_id, astral_sid, domain_start, domain_stop, astral_start, astral_stop, pdb_id, chain, overlap, session=None):
# Create a domain_to_astral ORM object and push it to the DB. Create a scoped session for push.
# Raise exception if push fails
    # If session not given, create scoped session for the DB push
    if not session:
        print "db_add:: No session provided. Creating scoped session"
        session = ScopedSession()
    
    if astral_start == None or astral_stop == None:
        dtoa_obj = DomainAstralOverlap(domain_id=domain_id, astral_id=astral_id, astral_sid=astral_sid, domain_start=domain_start, domain_stop=domain_stop, pdb_id=pdb_id, chain=chain, overlap=overlap)
    else:
        dtoa_obj = DomainAstralOverlap(domain_id=domain_id, astral_id=astral_id, astral_sid=astral_sid, domain_start=domain_start, domain_stop=domain_stop, astral_start=int(astral_start), astral_stop=int(astral_stop), pdb_id=pdb_id, chain=chain, overlap=overlap)
    session.add(dtoa_obj)
    try:
        session.flush()
    except IntegrityError:
        print "DomainToAstral {0} object already in database. Skipping".format(dtoa_obj)
        session.rollback()
    except Exception as e:
        print "Error in pushing DomainToAstral object {0} to database".format(dtoa_obj)
        raise
    return dtoa_obj


def main():
    # Call domain_to_astral
    experiment = 1171
    #storage = 'human_rna_astral_overlap.pkl'
    
    #domain_to_astral(experiment_id=experiment, store_file=storage, db_import=False)
    domain_to_astral(experiment_id=experiment, db_import=True)

if __name__ == "__main__":
    main()



