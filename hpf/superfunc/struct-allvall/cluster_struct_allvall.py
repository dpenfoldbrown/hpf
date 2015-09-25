#!/usr/bin/env python

# A script mimicing homolog_struct_allvall.py, but adapted to be superfast on the cluster
# No DB queries, instead holds everything in memory from file reads

import re
import time
import bisect
import tempfile
import cPickle as pickle
from hpf.utilities import consume
from hpf.processing import processor
from hpf.hddb.db import Session, ScopedSession, HomologStructAllVAll
from hpf.superfunc.homolog_struct_set import HomologComparisonProtein

# Global Setup

query_struct_pkl = '/scratch/superfunc/data/allhuman_largeblast_query_struct_dict.pkl'
num_pieces = 29
KEEP_BEST  = 100
RESULTS_DIR = '/tmp'

# Pickles storing the astral_index dict (astral_sid -> index int) and 
# astral all-v-all matrix (where list[i][j] is score of astrals with indices i and j)
astral_list_file   = '/scratch/superfunc/data/astral_mammoth.db.out'
astral_index_file  = '/scratch/superfunc/data/astral_index.pkl'
astral_matrix_file = '/scratch/superfunc/data/astral_matrix_empty.pkl'
astral_index  = None
astral_matrix = None


class HomologComparisonProteinPair():

    def __init__(self, source_hcprotein, target_hcprotein):
    # source and target are HomologComparisonProtein objects 
        self.source = source_hcprotein
        self.target = target_hcprotein
        self.similarity = self.similarity_score()

    def __lt__(self, other):
    # Reverse lessthan to facilitate ordered insertion from greatest to least.
        return self.similarity > other.similarity

    def __repr__(self, ):
        return "HomologCompProtPair: <Source {0}>, Target <{1}>, Similarity {2}".format(self.source, self.target, self.similarity)

    def similarity_score(self, ):
        total_score = 0.0
        total_score += self.pairwise_max_score(self.source.structure_list, self.target.structure_list)
        total_score += self.pairwise_max_score(self.target.structure_list, self.source.structure_list)

        return total_score / (self.source.num_domains + self.target.num_domains)

    def pairwise_max_score(self, source_list, target_list):
        score_sum = 0.0

        for source_struct in source_list:
            if source_struct == None:
                continue
            max_score = 0.0
            for target_struct in target_list:
                if target_struct == None:
                    continue
                try:
                    score = self._get_score(source_struct, target_struct)
                except KeyError as k:
                    print "Excepton: {0}. No score entry for astral {1} v astral {2}".format(k, source_struct, target_struct)
                    continue
                if score == None:
                    continue
                elif score > max_score:
                    max_score = score
            score_sum += max_score
        
        return score_sum

    def _get_score(self, src_struct, cmp_struct):
        if astral_matrix == None:
            raise("Astral all-v-all matrix has not been populated")
        if astral_index == None:
            raise("Astral index has not been populated")
        score = astral_matrix[astral_index[src_struct]][astral_index[cmp_struct]]
        if score == None:
            score = astral_matrix[astral_index[cmp_struct]][astral_index[src_struct]]
        return score


def output_file(outhandle, list):
# Takes a filehandle to output to (must be open/writable) and a list of
# HomologComparisonProteinPair objects to store in file. Form:
# source_queryid    target_queryid  similarity
    for pair in list:
        outhandle.write("{0}\t{1}\t{2}\n".format(pair.source.query_id, pair.target.query_id, pair.similarity))

def output_db(list):
# Takes a list of HomologComparisonProteinPair objects to store in the DB,
# table hpf.homolog_structure_allvall[...] 
    session = ScopedSession()
    for pair in list:
        hsa_obj = HomologStructAllVAll(source_id=pair.source.query_id, target_id=pair.target.query_id, score=pair.similarity)
        session.add(hsa_obj)
        try:
            session.flush()
        except IntegrityError:
            print "Object {0} already in database. Rolling back and skipping..".format(hsa_obj)
            session.rollback()
            continue
        except:
            print "Could not add protein comparison object {0} to DB".format(hsa_obj)
            raise

def output_pkl(outhandle, list):
# Takes a filehandle intended as a pickle to output to, and a list of HomologComparisonProteinPair 
# objects to store in pickle. Pickle entry: (source, target, similarity), a triple where source 
# and target are HomologComparisonProtein objects
    for pair in list:
        pickle.dump((pair.source, pair.target, pair.similarity), outhandle)


def do_allvall(query_ids):

    global astral_index, astral_matrix

    # Load astral index dict
    print "Loading astral index.."
    handle = open(astral_index_file)
    astral_index = pickle.load(handle)
    handle.close()

    # Load astral all-v-all score matrix columns to build up astral matrix
    #print "Loading and creating astral matrix.."
    #matrix_size = len(astral_index.keys())
    #astral_matrix = []
    #handle = open(astral_matrix_columns_file)
    #for i in range(matrix_size):
    #    column = pickle.load(handle)
    #    astral_matrix.append(column)
    #    print "{0}\r".format(i),
    #handle.close()

    # Load empty astral_matrix from pickle
    print "Loading empty astral matrix"
    handle = open(astral_matrix_file)
    astral_matrix = pickle.load(handle)
    handle.close()

    # Load astral all-v-all list and create astral matrix
    print "Creating astral_matrix from astral all-v-all list"
    start = time.time()
    handle = open(astral_list_file)
    astral_list_pattern = r"(?P<src_astral>[a-zA-Z0-9_]+)\s(?P<cmp_astral>[a-zA-Z0-9_]+)\s(?P<score>[0-9.]+)"
    for line in handle:
        found = re.match(astral_list_pattern, line)
        if not found:
            raise Exception("Line '{0}' doesn't match format".format(line))
        src_astral = found.group('src_astral')
        cmp_astral = found.group('cmp_astral')
        score = found.group('score')

        astral_matrix[astral_index[src_astral]][astral_index[cmp_astral]] = float(score)
        # DO NOT do reflection
        #astral_matrix[astral_index[cmp_astral]][astral_index[src_astral]] = float(score)
    handle.close()
    end = time.time()
    print "Astral matrix creation time (seconds): {0}".format(end-start)

    # Load query->HomologComparisonProtein dict
    print "Loading query structure dict.."
    handle = open(query_struct_pkl)
    query_dict = pickle.load(handle)
    handle.close()

    # For each query in task query set, do query v All struct comparison
    for source_query in query_ids:
        
        source_record = query_dict[source_query]
        print "Comparing query {0} (protein {1}) v. All".format(source_query, source_record.hit_protein)

        best_comparisons = []
        for cmp_record in query_dict.values():
            
            if source_record.query_id == cmp_record.query_id:
                continue
            pair_obj = HomologComparisonProteinPair(source_record, cmp_record)
            
            bisect.insort(best_comparisons, pair_obj)
            if len(best_comparisons) > KEEP_BEST:
                best_comparisons.pop()
            
        output_db(best_comparisons)
    
    print "Task complete"


def tasks(query_dict_file, num_pieces):
# Takes the filename storing the query -> HomologComparisonProtein dict
# and the number of pieces to divide the query IDs into to pass as a task
# to the allvall processing function
# A "task" represents a set of query IDs to process

    handle = open(query_dict_file)
    query_dict = pickle.load(handle)
    handle.close()

    query_ids = query_dict.keys()
    queries_per_task = len(query_ids) / num_pieces

    query_lists = []
    for i in range(0, len(query_ids), queries_per_task):
        query_lists.append(query_ids[i:i+queries_per_task])
    return query_lists

def main():
    print "Struct all-v-all main"
    pool = processor()
    pool.make_tasks(tasks, query_struct_pkl, num_pieces)
    consume(pool.run(do_allvall))

def test():
    print "Test run"
    handle = open(query_struct_pkl)
    query_dict = pickle.load(handle)
    handle.close()

    query_ids = query_dict.keys()
    do_allvall(query_ids[:1])

if __name__ == "__main__":
    #test()
    main()


