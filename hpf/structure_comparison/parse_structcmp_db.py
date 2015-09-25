"""
Functions to parse raw DB output (tab and line-delimited) for structure_mammoth and mammoth
DBs

@auth dpb
@date 1/07/2014
"""

import sys
from collections import namedtuple

StructPair = namedtuple("StructComp", ["source", "target"])

def parse_all():
    """Parses all files into same dict (see below for better code)"""
    #smfile00 = "/data/cafa/structMammothDB/structure_mammoth.txt"
    smfile01 = "/data/cafa/structMammothDB/structure_mammoth_01.txt"
    mfile00 = "/data/cafa/structMammothDB/mammoth.txt"

    struct_dict = {}
    #with open(smfile00) as handle:
    #    for line in handle:
    #        fields = line.split(",")
    #        struct_dict[StructPair(source=fields[1], target=fields[2])] = fields[10]
    #print "Completed parsing {0}".format(smfile00)
    with open(smfile01) as handle:
        for line in handle:
            fields = line.split(",")
            struct_dict[StructPair(source=fields[1], target=fields[2])] = fields[10]
    print "Completed parsing {0}".format(smfile01)
    with open(mfile00) as handle:
        for line in handle:
            fields = line.split(",")
            struct_dict[StructPair(source=fields[0], target=fields[1])] = fields[2]
    print "Completed parsing {0}".format(mfile00)
    return struct_dict

def parse_structure_mammoth(filename):
    """
    Parses structure_mammoth style tables:
        0 `id` int(11) NOT NULL AUTO_INCREMENT,
        1 `prediction_id` int(11) NOT NULL,
        2 `experiment_id` int(11) NOT NULL,
        3 `prediction_type` enum('astral','pdb','denovo','other') NOT NULL,
        4 `experiment_type` enum('astral','pdb','denovo','other') NOT NULL,
        5 `timestamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
        6 `ini_psi` double unsigned NOT NULL,
        7 `ini_rms` double DEFAULT NULL,
        8 `end_psi` double unsigned NOT NULL,
        9 `end_rms` double DEFAULT NULL,
        10 `zscore` double unsigned NOT NULL,
        11 `evalue` double unsigned NOT NULL,
        12 `version` int(11) NOT NULL,
    Returns dict: (source_id, target_id) => mammoth zscore
    (named tuple, float)
    """
    #count = 0
    struct_dict = {}
    with open(filename) as handle:
        for line in handle:
            #fields = line.rstrip().split(",")
            fields = line.split(",")
            #assert len(fields) == 13
            struct_dict[StructPair(source=fields[1], target=fields[2])] = fields[10]
            #count += 1
            #if count % 1000000 == 0:
            #    print "Done {0}".format(count)
    print "Completed parsing {0}".format(filename)
    print "Dict size: {0} MB".format(sys.getsizeof(struct_dict) / 1048576)
    return struct_dict

def parse_mammoth(filename):
    """
    Parses mammoth style tables:
        0 `p_structure_key` int(11) NOT NULL,
        1 `e_structure_key` int(11) NOT NULL,
        2 `zscore` double NOT NULL,
        3 `ini_psi` double unsigned NOT NULL,
        4 `end_psi` double unsigned NOT NULL,
        5 `score` double unsigned NOT NULL,
        6 `evalue` double unsigned NOT NULL,
        7 `ln_e` double unsigned NOT NULL,
        8 `num_p` int(11) unsigned NOT NULL,
        9 `num_e` int(11) unsigned NOT NULL,
        10 `nss` int(11) unsigned NOT NULL,
        11 `nsup` int(11) unsigned NOT NULL,
        12 `remark` text,
        13 `timestamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
        14 `ini_rms` double unsigned DEFAULT NULL,
        15 `end_rms` double unsigned DEFAULT NULL,
        16 `p_contact_order` double unsigned DEFAULT NULL,
        17 `e_contact_order` double unsigned DEFAULT NULL,
    Returns dict (source_id, target_id) => mammoth zscore
    (named tuple, float)
    """
    #count = 0
    struct_dict = {}
    with open(filename) as handle:
        for line in handle:
            #fields = line.rstrip().split(",")
            fields = line.split(",")
            #assert len(fields) == 18
            struct_dict[StructPair(source=fields[0], target=fields[1])] = fields[2]
            #count += 1
            #if count % 100000 == 0:
            #    print "Done {0}".format(count)
    print "Completed parsing {0}".format(filename)
    print "Dict size: {0} MB".format(sys.getsizeof(struct_dict) / 1048576)
    return struct_dict

if __name__ == "__main__":
    parse_all()
    #parse_structure_mammoth("/data/cafa/structMammothDB/structure_mammoth_01.txt")
    #parse_mammoth("/data/cafa/structMammothDB/mammoth.txt")

