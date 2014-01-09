"""
Functions for doing structure -> structure mammoth (via Mammoth commandline run)

@auth dpb 
@date 2/21/2013
"""

import os
from hpf.hddb.db import Session
session = Session()

mammoth_struct_type = {
    'astral'         : 'astral',
    'decoy'          : 'decoy',
    'homology_model' : 'other',
    'native'         : 'pdb',
    'native_groomed' : 'pdb',
    'native_pdbCompleted' : 'pdb',
    'native_processed'    : 'pdb',
    'pdbClean'            : 'pdb',
}

def sm_query(p_id, e_id, version):
    return _sm_query(p_id, e_id, version)


def _sm_query(p_id, e_id, version):
    """
    Queries both ways for a structure->structure mammoth:
    (p_id = X and e_id = Y) OR (p_id = Y and e_id = X)
    Only want to store one way in DB (avoid redundant storage), so have to double query
    (could put in all one query, but it looks like death)
    """
    from hpf.hddb.db import StructureMammoth, StructureMammoth01, StructureMammoth1186 
    
    #NOTE: Will now query all StructureMammoth tables (or, the base and 01, and all after) to
    #NOTE: attempt to get any SM record available. (Will not do org-specific tables here, bc 
    #NOTE: that was a bad idea in the first place). 
    #NOTE: Should in fact break tables into 50 or 100, where p_id % (mod) 50 denotes which table
    #NOTE: an entry should be stored in/queried for in. Write that into the db.py ORM code. 
    
    # StructureMammoth (original), normal and flipped
    score = session.query(StructureMammoth
        ).filter(
            (StructureMammoth.prediction_id==p_id) & 
            (StructureMammoth.experiment_id==e_id)
        ).filter(
            StructureMammoth.version==version
        ).first()
    if score:
        session.close()
        return score
    score = session.query(StructureMammoth
        ).filter(
            (StructureMammoth.prediction_id==e_id) & 
            (StructureMammoth.experiment_id==p_id)
        ).filter(
            StructureMammoth.version==version
        ).first()
    if score:
        session.close()
        return score
    
    # StructureMammoth01, normal and flipped
    score = session.query(StructureMammoth01
        ).filter(
            (StructureMammoth01.prediction_id==p_id) & 
            (StructureMammoth01.experiment_id==e_id)
        ).filter(
            StructureMammoth01.version==version
        ).first()
    if score:
        session.close()
        return score
    score = session.query(StructureMammoth01
        ).filter(
            (StructureMammoth01.prediction_id==e_id) & 
            (StructureMammoth01.experiment_id==p_id)
        ).filter(
            StructureMammoth01.version==version
        ).first()
    if score:
        session.close()
        return score
    
    return None

def _mammoth_query(p_id, e_id):
    """
    Queries both ways for a structure->structure mammoth in hpf.mammoth table
    NOTE: Currently returns None, as hpf.mammoth table stores all RMS values
    (ini_rms, end_rms) as NULL (none have value).
    """
    #from hpf.hddb.db import Session
    #session = Session()
    from hpf.hddb.db import Mammoth as M
    mammoth = session.query(M).filter((M.p_structure_key==p_id) & (M.e_structure_key==e_id)).first()
    if mammoth: 
        return mammoth
    mammoth = session.query(M).filter((M.p_structure_key==e_id) & (M.e_structure_key==p_id)).first()
    session.close()
    if mammoth:
        return mammoth
    return None

def write_struct_file(filename, structure):
    """Writes a PDB file of the given Structure DBO object to filename"""
    with open(filename, 'w') as handle:
        handle.write(structure.text)

def run_mammoth(prediction_id, experiment_id, ptype, etype):
    """Simply runs mammoth on two structures, returns results, doesn't store"""
    from hpf.hddb.db import Structure
    from hpf.pdb.mammoth import Mammoth, MammothCL

    print "RunMammoth"

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
    for file in files: 
        os.remove(os.path.join(work_dir,file))
    os.removedirs(work_dir)
    
    session.close()
    return mscore

def structure_mammoth(prediction_id, experiment_id, ptype, etype, base_dir=None, dbstore=True, version=1, 
        debug=False, cleanup=True, table_destination="default", fake_mammoth=False):
    """
    Takes the hpf DB structure IDs of two protein structures. Optionally checks for pre-computed
    mammoth in DB (if exists, returns). Fetches PDB files for both structures to local storage,
    runs mammoth on structures, and optionally stores mammoth results in DB (hpf, structure_mammoth tbl).
    
    Creates directory p<prediction_id>_e<experiment_id>_mammoth in given base_dir. Removes this directory
    if cleanup is True.

    version is optional if when specifying dbstore=True, you want to check for and add a specific version 
    (field in hpf.structure_mammoth) of structure->structure mammoth result. Defaults to 1
    """
    from hpf.pdb.mammoth import Mammoth, MammothAlignment, MammothCL
    from hpf.hddb.db import Session, Structure, Mammoth as MammothORM, StructureMammoth, StructureMammoth1186, StructureMammoth_class_lookup, push_to_db

    # Init session and check for pre-existing mammoth (optional)
    #session = Session()
    if dbstore:
        # Check StructureMammoth (hpf.structure_mammoth) for structure IDs
        previous = _sm_query(prediction_id, experiment_id, version)
        if previous:
            if debug:
                print "Previous version of Struct->Struct Mammoth exists: {0}".format(previous)
            return previous
        
        # Check Mammoth (hpf.mammoth) for structure IDs - if exists, transfer values, do not rerun Mammoth
        mammoth_dbo = _mammoth_query(prediction_id, experiment_id)
        if mammoth_dbo:
            if debug:
                print "Struct->Struct comparison exists in Mammoth table (hpf.mammoth): {0}".format(mammoth_dbo),
                print "Adding to StructureMammoth table and returning"
            SM = StructureMammoth_class_lookup[table_destination]
            sm_dbo = SM(
                prediction_id=mammoth_dbo.p_structure_key,
                experiment_id=mammoth_dbo.e_structure_key,
                prediction_type=ptype,
                experiment_type=etype,
                ini_psi=mammoth_dbo.ini_psi,
                ini_rms=mammoth_dbo.ini_rms,
                end_psi=mammoth_dbo.end_psi,
                end_rms=mammoth_dbo.end_rms,
                zscore=mammoth_dbo.zscore,
                evalue=mammoth_dbo.evalue,
                version=version)

            push_to_db(session, sm_dbo, raise_on_duplicate=False)
            if debug:
                print "Added: {0}".format(sm_dbo)
            return sm_dbo
    
    if fake_mammoth:
        #-------------------------------------------------------------------------------------
        # MANUAL BREAK: IF NO RESULTS IN THE DB, RETURN LARGE SCORE TO INDICATE MAMMOTH NECESSARY
        #print "Returning FALSE MAMMOTH"
        SM = StructureMammoth_class_lookup[table_destination]
        return SM(
           prediction_id=prediction_id,
           experiment_id=experiment_id,
           prediction_type=ptype,
           experiment_type=etype,
           ini_psi=0.0,
           ini_rms=0.0,
           end_psi=0.0,
           end_rms=0.0,
           zscore=100000.0,
           evalue=0.0,
           version=version)
        #-------------------------------------------------------------------------------------

    # Check input workdir
    if not base_dir:
        base_dir = os.getcwd()
    if not os.path.isdir(base_dir):
        raise Exception("Work dir '{0}' is not a valid directory".format(base_dir))
    
    # Get prediction and experiment struct DBOs
    p_structure = session.query(Structure).get(prediction_id)
    if not (p_structure):
        raise Exception("Prediction structure ({0}) not found in DB".format(prediction_id,))
    e_structure = session.query(Structure).get(experiment_id)
    if not (e_structure):
        raise Exception("Experiment structure ({0}) not found in DB".format(experiment_id,))

    # Set up working environment
    work_dir = os.path.join(base_dir, "p{0}_e{1}_mammoth".format(prediction_id, experiment_id))
    os.mkdir(work_dir)
    os.chdir(work_dir)
    prediction_file = "p{0}.pdb".format(prediction_id)
    experiment_file = "e{0}.pdb".format(experiment_id)
    write_struct_file(prediction_file, p_structure)
    write_struct_file(experiment_file, e_structure)

    # Run mammoth (via hpf.pdb.mammoth Mammoth), parse results
    mammoth_file = "p{0}_e{1}.mammoth".format(prediction_id, experiment_id)
    mcl = MammothCL(experiment_file, prediction_file, cwd=work_dir, output=mammoth_file)
    m = Mammoth(mcl, parse=True)
    mscore = m.run()

    # (Optional) store structure-structure mammoth in DB
    SM = StructureMammoth_class_lookup[table_destination]
    sm_dbo = SM(
        prediction_id=prediction_id,
        experiment_id=experiment_id,
        prediction_type=ptype,
        experiment_type=etype,
        ini_psi=mscore.ini_psi,
        ini_rms=mscore.ini_rms,
        end_psi=mscore.end_psi,
        end_rms=mscore.end_rms,
        zscore=mscore.zscore,
        evalue=mscore.evalue,
        version=version)

    if dbstore:
        push_to_db(session, sm_dbo, raise_on_duplicate=False)
    
    # Complete and (optional) cleanup
    print "Mammothing {0} against {1} complete (zscore {2})".format(prediction_id, experiment_id, mscore.zscore)
    os.chdir(base_dir)
    if cleanup:
        files = os.listdir(work_dir)
        for file in files: os.remove(os.path.join(work_dir,file))
        os.removedirs(work_dir)

    session.close()
    return sm_dbo

if __name__ == "__main__":
    import argparse, itertools as it

    parser = argparse.ArgumentParser(description="Structure->Structure Mammoth")
    
    parser.add_argument("-p", "--prediction_id", action="store", dest="p_id", 
            help="Predicted structure ID")
    parser.add_argument("-e", "--experiment_id", action="store", dest="e_id", 
            help="Experiment structure ID")
    parser.add_argument("--version", action="store", dest="version", default='1',
            help="Version of StructureMammoth in DB to check/store to")

    parser.add_argument("-f", "--structure_key_file", action="store", dest="structure_key_file", default="",
            help="File containing structure keys to be compared")
    
    parser.add_argument("--ptype", action="store", dest="p_type", required=True,
            help="Type of prediction structure (astral, pdb, denovo, other)")
    parser.add_argument("--etype", action="store", dest="e_type", required=True,
            help="Type of experiment structure (astral, pdb, denovo, other)")
    parser.add_argument("-d", "--dir", action="store", dest="dir", default=None,
            help="Base directory to put mammoth work files in")
    parser.add_argument("--mammoth_table", action="store", dest="table_destination", default="default",
            help="The StructureMammoth table in hpf database to store mammoth results, default=StructureMammoth, 1186=yeast experiment")
    
    parser.add_argument("--dbstore", action="store_true", dest="dbstore", default=False,
            help="Flag to store results in hpf.structure_mammoth table. Default False.")
    parser.add_argument("--cleanup", action="store_true", dest="cleanup", default=False,
            help="Flag to remove working files. Default False.")
    parser.add_argument("--verbose", action="store_true", dest="verbose", default=False,
            help="Flag to turn on debug messages. Default False.")

    args = parser.parse_args()

    if "" != args.structure_key_file:
        struct_keys = []
        handle = open(args.structure_key_file)
        for line in handle.readlines():
            try:
                int(line.rstrip())
                struct_keys.append(line.rstrip())
            except:
                continue



        for (p_id, e_id) in it.combinations(struct_keys,2):
            print "running: %s and %s" % (p_id, e_id,)

            structure_mammoth(p_id, 
                              e_id, 
                              args.p_type, 
                              args.e_type, 
                              base_dir=args.dir, 
                              dbstore=args.dbstore, 
                              version=args.version, 
                              cleanup=args.cleanup, 
                              debug=args.verbose,
                              table_destination=args.table_destination,
                             )    

    else:
        print "single mammoth run"
        structure_mammoth(args.p_id, 
                          args.e_id, 
                          args.p_type, 
                          args.e_type, 
                          base_dir=args.dir, 
                          dbstore=args.dbstore, 
                          version=args.version, 
                          cleanup=args.cleanup, 
                          debug=args.verbose,
                          table_destination=args.table_destination,
                         )    

