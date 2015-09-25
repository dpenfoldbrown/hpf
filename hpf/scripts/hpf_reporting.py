#!/usr/bin/env python

# Need to fill out and object-orientate. Avoid printing.
# A script to print data from proteins of a given experiment or set of sequences.
# A basic reporting script

from hpf.hddb.db import Session, Protein

def main():
    experiment = 1145

    session = Session()
    proteins = session.query(Protein).filter_by(experiment_key=experiment)
    #proteins = session.query(Protein).filter_by(experiment_key=experiment).limit(10)
    
    for protein in proteins:
        print ">hpf | {0} | {1} | {2} | {3}".format(protein.sequence_key, protein.id, protein.experiment_key, protein.sequence.sequence[:20]+"...")
        domains = get_old_domains(protein)
        if domains == [] or domains == None:
            continue
        count = 1
        for domain in domains:
            repr = "\tDomain {0}: {1}".format(count, domain.domain_type)
            if domain.sccs != None and domain.sccs.sccs != None:
                repr += "\t{0}".format(domain.sccs.sccs)
                if domain.sccs.scopDescription != None:
                    repr += " ({0})".format(domain.sccs.scopDescription.eng_desc)
            if domain.domain_type == 'psiblast' and domain.sccs != None:
                repr += "\t{0}:{1}".format(domain.sccs.pdbid, domain.sccs.chain)
            count += 1
            print repr
    print "Complete"

def get_old_domains(protein):
    if protein.all_domains == []:
        return protein.domains
    min_ginzu = protein.all_domains[0].ginzu_version
    for domain in protein.all_domains:
        if domain.ginzu_version < min_ginzu:
            min_ginzu = domain.ginzu_version

    oldest_domains = list()
    for domain in protein.all_domains:
        if domain.ginzu_version != min_ginzu:
            continue
        oldest_domains.append(domain)
    return oldest_domains

if __name__ == "__main__":
    main()
