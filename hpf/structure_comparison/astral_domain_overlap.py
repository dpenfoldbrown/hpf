#
# A script to compute the ratio of astral overlap between a ginzu domain of known structure
# and the corresponding PDB's component astral structures.
# Like this:
#
# Protein     |--------d1 (known: pdb A:50-100)------+------d2 (unknown)------|
#             |                                      | 
#             50                                     100
# PDB A   |...---------------------------------------------------------...|
# astrals |-----a1--|--------a2------|----------a3------|--------a4----...|
#
# a1, a2, and a3 are overlapping with domain d1 (aka PDB A:50-100)
#
# astral-domain overlap is the portion of astral (a1,a2,a3) that overlaps the domain region on
# matching PDB. EG:  a1-overlap = 6/9, a2-overlap = 1, a3-overlap = 16/18 
#
# dpb 2/21/2013
#


def astral_to_domain(experiment_id, threshold=0.5, dbstore=False):
    """
    Fetches all known-type domains for givein experiment ID. Computes astral overlap to
    domains, prints and optionally stores in hpf DB (table astral_domain_overlap)
    Will only store overlaps >= threshold parameter
    """
    from hpf.hddb.db import Session, push_to_db,  AstralDomainOverlap, Protein, Domain
    from hpf.structure_comparison.overlap import overlap
    from hpf.structure_comparison.astral_util import get_astrals, get_astral_startstop, parse_astral_chain

    # Create session and get all domains
    session = Session()
    domains = session.query(Domain).join(Protein).filter(Protein.experiment_key==experiment_id).filter(Domain.domain_type.in_(['psiblast','fold_recognition'])).all()
    print "Considering {0} domains to compute astral overlap for".format(len(domains))

    # For each domain, get representative astrals, calculate overlap, and store (optional)
    missing_astral = 0
    for domain in domains:
        domain_pdb_start = domain.region.parent_start
        domain_pdb_stop  = domain.region.parent_stop

        astrals = get_astrals(domain, session)
        if not astrals:
            #print "No astrals found for domain {0}".format(domain.id)
            missing_astral += 1
            continue

        for astral in astrals:
            try:
                (astral_start, astral_stop) = get_astral_startstop(astral)
                overlap_ratio = overlap(astral_start, astral_stop, domain_pdb_start, domain_pdb_stop)
            except ValueError:
                print "Negative overlap for domain {0} ({1}-{2}), Astral {3} (PDB {4}{5})".format(
                        domain.id, domain_pdb_start, domain_pdb_stop, astral.sid, astral.pdbid, astral.chain)
                print "Ignoring, moving to next astral.."
                continue
            except:
                print "Error calculating overlap for  domain {0} ({1}-{2}), Astral {3} (PBD {4}{5})".format(
                        domain.id, domain_pdb_start, domain_pdb_stop, astral.sid, astral.pdbid, astral.chain)
                raise

            if dbstore and overlap_ratio >= float(threshold):
                chain = parse_astral_chain(astral.chain)
                atod_dbo = AstralDomainOverlap(astral_id=astral.id, 
                                               astral_sid=astral.sid, 
                                               domain_id=domain.id, 
                                               astral_start=astral_start, 
                                               astral_stop=astral_stop, 
                                               domain_start=domain_pdb_start, 
                                               domain_stop=domain_pdb_stop, 
                                               pdb_id=astral.pdbid, 
                                               chain=chain, 
                                               overlap=overlap_ratio,
                                              )
                push_to_db(session, atod_dbo, exception_str="Error in pushing AstralDomainOverlap {0} to DB".format(atod_dbo), raise_on_duplicate=False)
            
            #print "Domain {0} ({1}-{2}), Astral {3} (PDB {4}{5}), Astral overlap {6}".format(domain.id, domain_pdb_start, domain_pdb_stop, astral.sid, astral.pdbid, astral.chain, overlap_ratio)

    print "Calculating astral to domain overlap for experiment {0} complete".format(experiment_id)
    print "{0} of {1} known-structure domains had no astral entries".format(missing_astral, len(domains))


if __name__ == "__main__":
    EXPERIMENT = 1202
    THRESHOLD  = 0.75
    DBSTORE    = True
    astral_to_domain(EXPERIMENT, threshold=THRESHOLD, dbstore=DBSTORE)

