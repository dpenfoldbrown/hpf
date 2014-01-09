#!/usr/bin/env python

# A script to drive the functionality of all-v-all structure comparison
# Does protein-protein structure comparison by astral structure coverage per
# protein domain (max average pairwise similarity over all domains).
# Works with a dictionary of for query_id -> HomologComparisonProtein object.
# Made to complement homolog_structure_set.py

import bisect
import tempfile
import cPickle as pickle
from hpf.utilities import consume
from hpf.processing import processor
from hpf.hddb.db import Session, ScopedSession, AstralComparison, HomologStructAllVAll
from hpf.superfunc.homolog_struct_set import HomologComparisonProtein


# Global setup
query_struct_pkl = '/home/dpb3/superfunc/allhuman_largeblast_query_struct_dict.pkl'
num_pieces = 500
KEEP_BEST = 100
RESULTS_DIR = '/tmp'


class HomologComparisonProteinPair():
# A class to represent a pair of HomologComparisonProteins, and to compute their
# similarity score. Mainly exists in order to define a __lt__ rich comparison method
# in order to be sortable by sim. score (and therefore inserted in-order w/bisect)

    def __init__(self, source_hcprotein, target_hcprotein):
    # source and target are HomologComparisonProtein objects
    # similarity is the similarity score between the two objects, based on the similarity_score function
        self.source = source_hcprotein
        self.target = target_hcprotein
        self.similarity = self.similarity_score()

    def __lt__(self, other):
        # This seems weird, but want to sort in greatest -> least order, so reverse less than
        return self.similarity > other.similarity

    def __repr__(self, ):
        return "HomologCompProtPair: <Source {0}>, Target <{1}>, Similarity {2}".format(self.source, self.target, self.similarity)
    
    def similarity_score(self, ):
        
        unweighted_score = 0.0
        
        # Get the sum of the max scores of structure similarity between source and target, and then
        # target and source (order is important, and need to do both ways)
        unweighted_score += self.pairwise_max_sum(self.source.structure_list, self.target.structure_list)
        unweighted_score += self.pairwise_max_sum(self.target.structure_list, self.source.structure_list)

        #DEBUG
        #print "Total unweighted score for {0} (prot: {4}) against {1} (prot: {5}) (with {2} total domains): {3}".format(self.source.query_id, \
        #        self.target.query_id, self.source.num_domains + self.target.num_domains, unweighted_score, self.source.hit_protein, self.target.hit_protein)
        
        return unweighted_score / (self.source.num_domains + self.target.num_domains)
    

    def pairwise_max_sum(self, source_list, target_list):
    # Returns the sum of the maximum pair score where a pair is a (top element, bottom element) pair
    # Implented to check the pdb.mammoth table ORM object, to get astral-v-astrall struct comp scores
    # If an element in a list is None, add nothing to total.
      
        #DEBUG
        #print "Source structs: ", source_list
        #print "Target structs: ", target_list

        session = ScopedSession()
        sum = 0.0
        
        # Compare each source structure to all target structures
        for source_struct in source_list:
            if source_struct == None:
                continue
            max_score = 0.0
            for target_struct in target_list:
                if target_struct == None:
                    continue
                struct_comp_obj = session.query(AstralComparison).filter_by(prediction=source_struct, experiment=target_struct).first()
                if not struct_comp_obj:
                    struct_comp_obj = session.query(AstralComparison).filter_by(prediction=target_struct, experiment=source_struct).first()
                    if not struct_comp_obj:
                        #DEBUG
                        #print "pms:: could not fetch comparison from the DB for target: {0}".format(target_struct)
                        continue
                #DEBUG
                #print "pms:: struct comp object found for target: {0}, score {1}, cur max: {2}".format(target_struct, float(struct_comp_obj.zscore), max_score)
                
                if  float(struct_comp_obj.zscore) > max_score:
                    max_score = float(struct_comp_obj.zscore)
            
            #DEBUG
            #print "Max score: {0} for source '{1}' against targets:".format(max_score, source_struct), target_list
            
            sum += max_score
        return sum


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
# Takes a list of query IDs. Computes structural similarity between each given ID
# and every record in the global results file containing query_id->HomologComparisonProtein
# dict.
# References the global query_struct_pkl to open dict of full query results

    #DEBUG
    print "Query IDs for all-v-all", query_ids

    # Create an outfile via tempfile, and open it for writing
    (fd, outfile) = tempfile.mkstemp(suffix=".out", prefix="human_allvall_", dir=RESULTS_DIR)
    outfile_handle = open(outfile, 'w')

    # Load dict of query->HomologComparisonProtein objects
    handle = open(query_struct_pkl)
    query_dict = pickle.load(handle)
    handle.close()

    # For each query, do query v all queries analysis
    for source_query in query_ids:
        source_record = query_dict[source_query]
        #DEBUG
        print "Comparing query {0} (protein {1}) v. All".format(source_query, source_record.hit_protein)
        
        best_similarities = []
        for cmp_record in query_dict.values():
            # Skip comparing against self
            if source_record.query_id == cmp_record.query_id:
                continue
            
            #DEBUG
            #print ""
            #print "Finding similarity of source {0} v target {1}".format(source_query, cmp_record.query_id)
            
            pair_obj = HomologComparisonProteinPair(source_record, cmp_record)
            bisect.insort(best_similarities, pair_obj)
            if len(best_similarities) > KEEP_BEST:
                best_similarities = best_similarities[:KEEP_BEST]

        # Output results list of best similarity scores between query and all
        output_file(outfile_handle, best_similarities)
        #output_db(best_similarities)
        
        #DEBUG
        #print "Query {0} best {1} similar queries/structs: ".format(source_query, KEEP_BEST)
        #for sim in best_similarities:
        #    print "\t{0}".format(sim)


    outfile_handle.close()
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
    global query_struct_pkl, num_pieces

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
    test()
    #main()


