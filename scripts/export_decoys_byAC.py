##
## A quick script to dump de novo structure filse for all HUMAN (experiments 804, 1176)
## protein domains.
## Identified via filename: <uniprot ac>_<protein#>_<ginzu#>_<domain#>_<domain_start>-<domain_stop>.pdb
##
## dpb 2/04/2013
##

from hpf.hddb.db import Session, Protein

# Globals
EXPERIMENTS = [1176, 804, ]

# Init DB
session = Session()

# Fetch all human proteins
print "Obtaining all proteins from Experiments {0}".format(EXPERIMENTS)
proteins = session.query(Protein).filter(Protein.experiment_key.in_(EXPERIMENTS)).all()
print "{0} proteins found".format(len(proteins))

# Consider all proteins. Export struct files for those with structure
written = 0
for p in proteins:
    # Skip entries with no AC membership
    if not p.ac:
        continue

    # For all domains of the protein, new and old (versions of ginzu)
    for d in p.all_domains:
        
        # If d has been chosen as foldable
        if d.foldable:
            
            # If d has been clustered and stored
            if d.foldable.convergence:
                
                # Get convergence record for cluster run on most decoys
                convergence_record = sorted(d.foldable.convergence, key=lambda r: r.total_decoys, reverse=True)[0]

                # Check for cluster centers
                if convergence_record.cluster_centers:
                    # Grab the first cluster center struct (first is biggest - already sorted)
                    center_struct = convergence_record.cluster_centers[0]

                    # Output it
                    outfile = "{0}_{1}_{2}_{3}_{4}-{5}.pdb".format(p.ac.ac, p.id, d.ginzu_version, d.domain_nr, d.region.start, d.region.stop)
                    handle = open(outfile, 'w')
                    handle.write(center_struct.text)
                    handle.close()
                    print "Wrote structure file {0}".format(outfile)
                    written += 1

print "Complete. Wrote {0} structure files".format(written)


            


