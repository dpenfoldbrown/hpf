"""
All functionality, in one big-ass function, globals, whole shebangerang, of protein-
protein structure similarity matrix code. Because it don't get no better.

@auth dpb
@date 1/09/2014

NOTE: This is the "IN MEMORY" version
NOTE: Using a global protein structure set cacheing dict (no re-building repr. struct sets)
"""

## CONFIG

EXPERIMENT = 1198
MTX_OUTFILE = "/home/dpb/hpf/hpf/structure_comparison/out/1198_mycge_mtx.out"
KEY_OUTFILE = "/home/dpb/hpf/hpf/structure_comparison/out/1198_mycge_key.out"
TOMAMMOTH_FILE = "/home/dpb/hpf/hpf/structure_comparison/out/1198_mycge_2mammoth.out"


## IMPORTS
import os
from collections import namedtuple

from hpf.pdb.mammoth import Mammoth, MammothCL
from hpf.hddb.db import Session, Protein, Structure

from hpf.structure_comparison.domain_structure_representation import structure_representation


## GLOBALS
session = Session()
StructPair = namedtuple("StructComp", ["source", "target"])

global_struct_mammoth_db = {}
global_protein_struct_sets = {}
global_org_structs = set()


## MAIN FUNCTIONALITY
def create_structure_matrix(experiment, mtxfile, keyfile):
    """All functionality for getting data and creating matrix outfile"""

    # Get list of all protein DBOs
    proteins = session.query(Protein).filter_by(experiment_key=experiment).all()
    if not proteins:
        raise Exception("Error: no proteins for experiment ID {0}".format(experiment))
    print "{0} proteins found in experiment {1}".format(len(proteins), experiment)

    # Loop over prots. Write in order to keyfile; build org struct set; build prot structset dict
    keyhandle = open(keyfile, 'w')
    keyhandle.write("UNIPROT\tHPF_SEQ\tHPF_PROT\n")
    
    for p in proteins:
        # Write entry to keyfile (uniprot ac, hpf seq id, hpf prot id)
        uniprotac = p.ac.ac if p.ac else " "*5
        keyhandle.write("{0}\t{1}\t{2}\n".format(uniprotac, p.sequence_key, p.id))

        # Get protein struct set, add it to global protein structset dict
        struct_set = get_structset(p)
        global_protein_struct_sets[p.id] = struct_set

        # Add structs in struct set to organism-level struct set
        # NOTE: Must loop here because struct_set can have nested lists
        for sv in struct_set:
            if is_known(sv): 
                global_org_structs.add(str(sv))
            else:
                global_org_structs.update(map(str, sv))
    
    keyhandle.close()
    print "{0} structures in Organism Structure set".format(len(global_org_structs))

    # Parse and store StructureMammoth DB files (store iff both structs in org struct set)
    print "IMPORTANT! Parsing StrucuteMammoth DB files. May take a while..."
    parse_structure_mammoth_db()
    print "Parsing StructureMammoth DB files complete"

    # Iterate symmetrically over list: for each pair i,j, compute i,j similarity,
    # store in matrix row for i, write row i to file, del row, move on to next row
    mtxhandle = open(mtxfile, 'w')
    num_proteins = len(proteins)
    i = 0
    while i < num_proteins:
        print "Processing row {0} of {1}".format(i+1, num_proteins)
        row = [0]*num_proteins
        
        # Get the outer-values only once per row.
        psource = proteins[i]
        source_structset = global_protein_struct_sets[psource.id]
        
        j = i + 1
        while j < num_proteins:
            ptarget = proteins[j]
            target_structset = global_protein_struct_sets[ptarget.id]

            # Pairwise max average
            forward_scores = unidirectional_pma(source_structset, target_structset, i, j)
            backward_scores = unidirectional_pma(target_structset, source_structset, i, j)
            scores = forward_scores + backward_scores 
            
            if len(scores) < 1:
                # Set row[j] to 0 and move to next protein
                row[j] = 0.0
                j += 1
                continue
            
            # Sum and average to get bi-directional PMA score, set row[j]
            row[j] = sum(scores) / len(scores)
            j += 1
            
        # Write row to matrix outfile
        for v in row:
            mtxhandle.write("{0}\t".format(v))
        mtxhandle.write("\n")
        mtxhandle.flush()
        
        i += 1
        print "Completed row {0}".format(i)

    mtxhandle.close()
    global_tomammoth_handle.close()
    print "Writing Matrix and Keyfile for experiment {0} COMPLETE".format(experiment)


## SUPPORT FUNCTIONS

def parse_structure_mammoth_db():
    """Parse in all structure mammoth db files, store in global_struct_mammoth_db dict"""
    smfile00 = "/data/cafa/structMammothDB/structure_mammoth.txt"
    smfile01 = "/data/cafa/structMammothDB/structure_mammoth_01.txt"
    mfile00 = "/data/cafa/structMammothDB/mammoth.txt"

    struct_dict = {}
    with open(smfile00) as handle:
        for line in handle:
            fields = line.split(",")
            if fields[1] in global_org_structs and fields[2] in global_org_structs:
                #print "Adding Structure Mammoth to global dict"
                global_struct_mammoth_db[StructPair(source=int(fields[1]), target=int(fields[2]))] = float(fields[10])
    print "Completed parsing {0}".format(smfile00)
    
    with open(smfile01) as handle:
        for line in handle:
            fields = line.split(",")
            if fields[1] in global_org_structs and fields[2] in global_org_structs:
                #print "Adding Structure Mammoth to global dict"
                global_struct_mammoth_db[StructPair(source=int(fields[1]), target=int(fields[2]))] = float(fields[10])
    print "Completed parsing {0}".format(smfile01)

    with open(mfile00) as handle:
        for line in handle:
            fields = line.split(",")
            if fields[0] in global_org_structs and fields[1] in global_org_structs:
                #print "Adding Structure Mammoth to global dict"
                global_struct_mammoth_db[StructPair(source=int(fields[0]), target=int(fields[1]))] = float(fields[2])
    print "Completed parsing {0}".format(mfile00)


def unidirectional_pma(source_structset, target_structset, protein_i, protein_j):
    """All-v-all struct sim. on given sets. Returns list of max scores for each pair"""
    max_pair_scores = []
    for sstruct in source_structset:
        source_scores = []
        for tstruct in target_structset:
            source_score = 0
            if is_known(sstruct):
                if is_known(tstruct):       # BOTH KNOWN
                    source_score = mem_structure_similarity(sstruct, tstruct, protein_i, protein_j)
                else:                       # Source KNOWN, Target UNKNOWN
                    max_pair_score = 0
                    for tstruct_decoy in tstruct:
                        score = mem_structure_similarity(sstruct, tstruct_decoy, protein_i, protein_j)
                        if score > max_pair_score:
                            max_pair_score = score
                    source_score = max_pair_score
            else:
                if is_known(tstruct):       # Source UNKNOWN, Target KNOWN
                    max_pair_score = 0
                    for sstruct_decoy in sstruct:
                        score = mem_structure_similarity(sstruct_decoy, tstruct, protein_i, protein_j)
                        if score > max_pair_score:
                            max_pair_score = score
                    source_score = max_pair_score
                else:                       # BOTH UNKNOWN
                    max_pair_score = 0
                    for sstruct_decoy in sstruct:
                        for tstruct_decoy in tstruct:
                            score = mem_structure_similarity(sstruct_decoy, tstruct_decoy, protein_i, protein_j)
                            if score > max_pair_score:
                                max_pair_score = score
                    source_score = max_pair_score
            source_scores.append(source_score)
        # Take the max of all the scores for this source structure
        if len(source_scores) > 0:
            max_source_score = max(source_scores)
            max_pair_scores.append(max_source_score)     
    return max_pair_scores

def is_known(structval):
    """Return True if structval is known type (single key), False if unknown type (list of keys)"""
    return not isinstance(structval, list)

def mem_structure_similarity(source, target, protein_i, protein_j):
    """In-memory version of structure similarity - get from dict or write to file"""
    sp = StructPair(source=source, target=target)
    if sp in global_struct_mammoth_db:
        #print "Found source {0}, target {1} in memory".format(source, target)
        return global_struct_mammoth_db[sp]

    sp = StructPair(source=target, target=source)
    if sp in global_struct_mammoth_db:
        #print "Found source {0}, target {1} in memory".format(target, source)
        return global_struct_mammoth_db[sp]

    # WRITE TOMAMMOTH FILE
    #print "Structs {0}, {1} not in memory, writing to file and skipping".format(source, target)
    global_tomammoth_handle.write("{0},{1},{2},{3}\n".format(protein_i, protein_j, source, target))
    global_tomammoth_handle.flush()
    return 0.0

    # RUN MAMMOTH (uncomment, and comment WRITE TOMAMMOTH FILE section above)
    #print "Structs {0}, {1} not in memory, mammothing".format(source, target)
    #mammoth_result = run_mammoth(source, target)
    #if not mammoth_result:
    #    raise Exception("RunMammoth: no score for structures {0}, {0}".format(source, target))
    #return float(mammoth_result.zscore)

def run_mammoth(prediction_id, experiment_id):
    """Runs mammoth on two structures, returns results object, doesn't store"""

    p_structure = session.query(Structure).get(prediction_id)
    assert p_structure
    e_structure = session.query(Structure).get(experiment_id)
    assert e_structure

    base_dir = os.getcwd()
    work_dir = os.path.join(base_dir, "p{0}_e{1}_mammoth".format(prediction_id, experiment_id))
    os.mkdir(work_dir)
    os.chdir(work_dir)
    prediction_file = "p{0}.pdb".format(prediction_id)
    experiment_file = "e{0}.pdb".format(experiment_id)
    write_struct_file(prediction_file, p_structure)
    write_struct_file(experiment_file, e_structure)

    mammoth_file = "p{0}_e{1}.mammoth".format(prediction_id, experiment_id)
    mcl = MammothCL(experiment_file, prediction_file, cwd=work_dir, output=mammoth_file)
    m = Mammoth(mcl, parse=True)
    mscore = m.run()

    os.chdir(base_dir)
    files = os.listdir(work_dir)
    for f in files:
        os.remove(os.path.join(work_dir, f))
    os.removedirs(work_dir)

    session.close()
    return mscore

def write_struct_file(filename, structure):
    """Writes a PDB file of the given Structure DBO object to filename"""
    with open(filename, 'w') as handle: 
        handle.write(structure.text)

def get_structset(protein):
    """Get representative structure set for protein. Returns list of ints"""
    struct_set = []
    for domain in protein.domains:
        domain_structure_set = structure_representation(domain)
        if not domain_structure_set:
            continue
        elif domain.known_type:
            struct_set += domain_structure_set
        else:
            struct_set.append(domain_structure_set)
    return struct_set


## MAIN and EXECUTION

if __name__ == "__main__":
    
    global_tomammoth_handle = open(TOMAMMOTH_FILE, 'w')
    global_tomammoth_handle.write("PROTEIN_I,PROTEIN_J,SOURCEID,TARGETID\n")

    create_structure_matrix(EXPERIMENT, MTX_OUTFILE, KEY_OUTFILE)




