
from hpf.hddb.db import Session, Protein
from hpf.structure_comparison.domain_structure_representation import *

session = Session()

# Get all proteins for experiment
experiment = 1186   # New yeast
proteins = session.query(Protein).filter_by(experiment_key=experiment).all()

# For all proteins, check domain struct. rep.
print "Experiment {0}, Proteins {1}".format(experiment, len(proteins))
empty = 0
for p in proteins:
    if p.id == 1523333:
        print "******"
        break 
    print ".",
    structure_keys = []
    for d in p.domains:
        domain_struct_set = structure_representation(d)
        if not domain_struct_set:
            continue
        elif d.known_type:
            structure_keys += domain_struct_set
        else:
            structure_keys.append(domain_struct_set)
    if structure_keys:
        print "Protein {0}, Structures {1}".format(p, structure_keys)
    else:
        empty += 1

print "{0} of {1} proteins have no domain structures".format(empty, len(proteins))
