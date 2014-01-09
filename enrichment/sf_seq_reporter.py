#!/usr/bin/env python

# A script to take a list of sequence keys (from file) or an experiment ID
# and a target SCOP superfamily (a.1.1), and report all the protein sequence IDs
# that have domains that have target structures of target SF

import sqlalchemy
from hpf.hddb.db import url
from sqlalchemy.sql import func, and_

engine = sqlalchemy.create_engine(url+'pwinters')
metadata = sqlalchemy.MetaData(engine)
pdomainsccs = sqlalchemy.Table('p_domain_sccs',metadata, autoload=True)

hpf_engine = sqlalchemy.create_engine(url+'hpf')
hpf_metadata = sqlalchemy.MetaData(hpf_engine)
scop_des = sqlalchemy.Table('scop_des',hpf_metadata, autoload=True)
seq_genesym = sqlalchemy.Table('sequenceGeneName', hpf_metadata, autoload=True)


def main():
    #superfams = ['d.68.6',]
    # All enriched SFs (for ~770 human rnab proteins against all-human background) with corrected pval >= 0.05
    superfams = ['d.58.7', 'd.51.1', 'i.1.1', 'a.21.1', 'b.40.4', 'd.50.1', 'd.64.1', 'b.38.1', 'c.37.1', 'b.131.1', 'a.4.5', 'd.68.6', 'c.120.1', 'a.289.1', 'd.79.3', 'd.43.1', 'c.51.1', 'a.118.19', 'd.265.1']
    
    seq_file = "/Users/dpb/bonneau-dev/hpf/trunk/src/projects/berliner_mrna/sequence_key_allquant_nonredundant.txt"
    seqs = parse_seqs(seq_file)
    sf_seqs_dict = get_superfamily_predictions_seqkey(seqs)
    
    # Human RNA: experiments = [1177,]
    #experiments = [804,]
    #sf_seqs_dict = get_superfamily_predictions(experiments)
    
    for target_sf in superfams:
        sf_seqs = sf_seqs_dict[target_sf]
        print "Superfamily {0} {1}:".format(target_sf, scopDesc(target_sf))
        print "Number of sequences: {0}".format(len(sf_seqs))
        for seq, domain_type in sf_seqs:
            print "{0}\t{1}\t{2}".format(seq, domain_type[:8], get_genesym(seq))
        
    print "Complete"    

def get_genesym(seq_key):
# Attempt to pull gene symbol from DB
    query = sqlalchemy.select([seq_genesym.c.symbol], seq_genesym.c.sequence_key == int(seq_key))
    sym = query.execute()
    sym_str = sym.fetchone()
    if not sym_str:
        return "None found"
    return sym_str[0]

def parse_seqs(file):
    seqs = list()
    handle = open(file)
    for line in handle:
        seqs.append(line.rstrip())
    handle.close()
    return seqs
    
def scopDesc(sccs):
    squery = sqlalchemy.select([scop_des.c.eng_desc], scop_des.c.sccs == sccs) 
    desc = squery.execute()
    return desc.fetchone()

def family2superfamily(sccs):
    return '.'.join(sccs.strip().split('.')[:-1])

def family2fold(sccs):
    return '.'.join(sccs.strip().split('.')[:-2])

def get_superfamily_predictions_seqkey(sequence_keys,prob_min=0.8):
    sf_dict = dict()
    squery = sqlalchemy.select([pdomainsccs.c.sccs,pdomainsccs.c.domain_type,pdomainsccs.c.parent_sequence_key], sqlalchemy.and_(pdomainsccs.c.sccs != "NULL", pdomainsccs.c.parent_sequence_key.in_(sequence_keys), pdomainsccs.c.confidence >= prob_min) )
    preds = squery.execute()
    for pred in preds:
        #print "sequence_key prediction %s" % pred
        for sccs in pred[0].split(','):
            try:
                sf_dict[ family2superfamily(sccs) ].append((pred[2], pred[1]))
            except KeyError:
                sf_dict[ family2superfamily(sccs) ] = [(pred[2], pred[1])]
    return sf_dict

def get_superfamily_predictions(experiments,prob_min=0.8):
    sf_dict = dict()
    squery = sqlalchemy.select([pdomainsccs.c.sccs,pdomainsccs.c.domain_type,pdomainsccs.c.parent_sequence_key], sqlalchemy.and_(pdomainsccs.c.sccs != "NULL", pdomainsccs.c.id.in_(experiments), pdomainsccs.c.confidence >= prob_min) )
    preds = squery.execute()
    for pred in preds:
        #print "sequence_key prediction %s" % pred
        for sccs in pred[0].split(','):
            try:
                sf_dict[ family2superfamily(sccs) ].append((pred[2], pred[1]))
            except KeyError:
                sf_dict[ family2superfamily(sccs) ] = [(pred[2], pred[1])]
    return sf_dict

if __name__ == "__main__":
    main()