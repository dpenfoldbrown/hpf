#
# Quick script to translate all entries in the astal_mammoth table into structure_mammoth
# entries (identified by structure IDs, not astral sids)
#
# dpb 2/21/2013
#

# TOO BIG TO DO IN ONE RUN (because sqlalchemy actually sucks... at iterated queries with MySQL
# (they don't work... still loads entire query result into memory).
# Do in pieces. There are:
# 142,537,733 
# entries. Split, and do in manageable-sized pieces

from hpf.hddb.db import Session, push_to_db, Astral, AstralComparison, StructureMammoth


def _sm_query(p_id, e_id, version, session):
    """
    Queries both ways for a structure->structure mammoth:
    (p_id = X and e_id = Y) OR (p_id = Y and e_id = X)
    Only want to store one way in DB (avoid redundant storage), so have to double query
    (could put in all one query, but it looks like death)
    """
    from hpf.hddb.db import StructureMammoth as SM
    
    # Normal
    score = session.query(SM).filter((SM.prediction_id==p_id) & (SM.experiment_id==e_id)).filter(SM.version==version).first()
    if score:
        return score
    # Flipped
    score = session.query(SM).filter((SM.prediction_id==e_id) & (SM.experiment_id==p_id)).filter(SM.version==version).first()
    session.close()
    if score:
        return score
    return None


ENTRIES = 142537733
PIECE_SIZE = 1000
RUNS = (ENTRIES / PIECE_SIZE) + 1


VERSION = 1
print "Translating astral_mammoth entries to structure_mammoth (VERSION {0})".format(VERSION)

# Start session
session = Session()

# For each entry, get struct ids, check if in structure_mammoth, add if not (use iterator query)
added = 0

START = 56662001
STOP  = START + PIECE_SIZE 
for i in range(RUNS):

    print "Translating records: {0}\t{1}".format(START, STOP)
    
    for ac in session.query(AstralComparison).filter(AstralComparison.cmp_id >= START).filter(AstralComparison.cmp_id < STOP).all():
     
        # Ignore if prediction (source) and comparison (experiment) astrals not present
        if not (ac.source_astral and ac.comparison_astral):
            continue
    
        # Check for pre-exisiting entry (either way, p->e or e->p) in StructureMammoth
        previous = _sm_query(ac.source_astral.structure_key, ac.comparison_astral.structure_key, VERSION, session)
        if previous:
            print "Mammoth structure {0} -> {1} already exists in DB".format(ac.source_astral.structure_key, ac.comparison_astral.structure_key)
            continue

        # Make StructureMammoth object and push to DB
        sm_dbo = StructureMammoth(prediction_id=ac.source_astral.structure_key,
                                  experiment_id=ac.comparison_astral.structure_key,
                                  prediction_type='astral',
                                  experiment_type='astral',
                                  ini_psi=ac.ini_psi,
                                  ini_rms=ac.ini_rms,
                                  end_psi=ac.end_psi,
                                  end_rms=ac.end_rms,
                                  zscore=ac.zscore,
                                  evalue=ac.evalue,
                                  version=VERSION
                                 )
        push_to_db(session, sm_dbo, raise_on_duplicate=False)
        added += 1

    START += PIECE_SIZE
    STOP  += PIECE_SIZE
    with open("astralmammoth2structuremammoth_workingon.txt", 'w') as handle:
        handle.write("AstralMammoth to StructureMammoth conversion: Working on piece {0} to {1}\n".format(START, STOP))

print "Translated {0} entries from astral_mammoth to structure_mammoth".format(added)
print "Complete"
    

