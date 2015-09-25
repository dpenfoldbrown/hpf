#!/usr/bin/env python
# Analyse the since solved stuffs from the DB
# A bunch of goofy bining and sorting to look at it the way I want

import sys
import MySQLdb
import getopt
from collections import defaultdict

def sf(sccs):
    split = sccs.split(".")[:3]
    if len(split) < 3:
        raise Exception("SCCS isn't long enough",sccs)
    return ".".join(split)

def prob(num, denom, digits=1):
    return str(round(div(num,denom)*float(100.0),digits))+"%"

def div(num, denom):
    if float(denom) == 0:
        return 0
    else:
        return float(num)/float(denom)
    
def top_scores(scores, N=None):
    """Grab the top N scores per unique id"""
    all = []
    
    # Bin by id
    id_scores = defaultdict(lambda: [])
    for p,m,s,i,c in scores:
        id_scores[i].append((p,m,s,i,c))

    # Pick N-best per ID    
    for id in id_scores:
        scores = id_scores[id]
        scores.sort(reverse=True)
        if N == None:
            N=len(scores)
        all+=scores[:N]
    
    all.sort(reverse=True)
    return all

def threshold(list, cutoff):
    """Threshold the list."""
    all = []
    for p,m,s,i,c in list:
        if p>= cutoff:
            all.append((p,m,s,i,c))
    return all

def correct(list):
    """Return all the scores that were correct"""
    cor = []
    for p,m,s,i,c in list:
        if sf(m)==sf(s):
            cor.append((p,m,s,i,c))
    return cor

def unique_ids(scores):
    """Return a set of the unique ids in the list"""
    return set([i for p,m,s,i,c in scores])

def unique_clusters(scores):
    """Return a set of the unique clusters in the list"""
    return set([c for p,m,s,i,c in scores])

def unique_sf(scores):
    return set([sf(s) for p,m,s,i,c in scores])

def subset_clusters(scores, *clusters):
    """Return all scores belonging to given clusters"""
    all = []
    for p,m,s,i,c in scores:
        if c in clusters:
            all.append((p,m,s,i,c))
    return all

def subset_sf(scores, *superfamilies):
    """Return all scores belonging to given superfamily"""
    #print superfamilies,len(scores)
    all = []
    for p,m,s,i,c in scores:
        if sf(s) in superfamilies:
            all.append((p,m,s,i,c))
    return all

def subset_class(scores, *scop_class):
    """Return all scores belonging to given scop class"""
    all = []
    for p,m,s,i,c in scores:
        if s[0] in scop_class:
            all.append((p,m,s,i,c))
    return all

def cutoffs(list):
    all = []
    for x in list:
        cutoff = div(x,10)
        all.append(cutoff)
    return all
                
if __name__=="__main__":
    cdhit_name = sys.argv[1]
    db = MySQLdb.connect(host="mcpeepants.bio.nyu.edu", user="patrick", passwd="patrick_nyu", db="since_solved")
    cursor = db.cursor()
    
    for table,column,index in [("mcmData","experiment_sccs","sequence_key"),("mcm_data_adjusted","sccs","domain_sequence_key")]:
        # A bunch of bins and thresholds for mcm scores
        query = "select s.cluster_id, s.sccs, c.identifier, m.probability, m.%s from since_solved s join cdhit_clstr c join hddb.%s m on s.cluster_id=c.cluster_id and c.identifier=m.%s where c.db='%s' and s.db='%s'" % (column, table, index, cdhit_name, cdhit_name)
        domain_query = "select s.cluster_id, s.sccs, c.identifier, m.probability, m.%s from since_solved s join cdhit_clstr c join hddb.%s m join hddb.domain d on s.cluster_id=c.cluster_id and c.identifier=m.%s and m.%s=d.domain_sequence_key where c.db='%s' and s.db='%s' and d.domain_type in ('psiblast','fold_recognition','pfam')" % (column, table, index, index, cdhit_name, cdhit_name)
        print query
        cursor.execute(query)
        # [cluster][hpf_id] yields list of scores
        all_scores = []
        #clusters = defaultdict(lambda: defaultdict(lambda: []))
        for cluster_id, sccs, hpf_id, probability, mcm_sccs in cursor.fetchall():
            # Group by cluster and then by hpf_id
            all_scores.append((probability, mcm_sccs, sccs, hpf_id, cluster_id))
    
        print "==",table,"=="
        
        if table=="mcmData":
            l = [8,9]
        else:
            l = [8]
        for cutoff in cutoffs(l):
            if table=="mcmData":
                l = [(1,"best1"),(None,"all5")]
            else:
                l = [(None,"all5")]
            for consider,name in l:
                print "Probability>=",cutoff,"using",name
                scores = top_scores(threshold(all_scores,cutoff),consider)
                print "\t",len(unique_clusters(scores)),"clusters with",len(unique_ids(scores)),"total sequences"
                print "\t",len(unique_sf(scores)),"unique superfamilies in threshold"
                correct_scores = correct(scores)
                
                count_clusters = len(unique_clusters(scores))
                count_correct_clusters = len(unique_clusters(correct_scores))
                print "\t",count_correct_clusters,"/",count_clusters," clusters correct:",prob(count_correct_clusters,count_clusters)
                print "\t",len(correct_scores),"/",len(scores)," mcm scores correct:",prob(len(correct_scores),len(scores))
                print "\t",len(unique_sf(correct_scores)),"unique superfamilies correct"
                for scop_class in ["a","b","c","d","e","f","g"]:
                    class_scores = subset_class(scores,scop_class)
                    class_scores_correct = subset_class(correct_scores,scop_class)
                    if len(class_scores) > 0:
                        print "\t\t",scop_class,":",len(unique_clusters(class_scores_correct)),"/",len(unique_clusters(class_scores)),"clusters:",prob(len(unique_clusters(class_scores_correct)),len(unique_clusters(class_scores)))
                        print "\t\t   ",len(unique_ids(class_scores_correct)),"/",len(unique_ids(class_scores)),"sequences:",prob(len(unique_ids(class_scores_correct)),len(unique_ids(class_scores)))
        
        print ""
#        
#        for superfamily in sorted(unique_sf(all_scores)):
#            scores = subset_sf(all_scores,superfamily)
#            one = top_scores(scores, 1)
#            all = top_scores(scores)
#            print "Superfamily",superfamily,"has",len(unique_clusters(scores)),"clusters with",len(unique_ids(scores)),"sequences"
#            for cutoff in cutoffs([8,9]):
#                if len(threshold(one, cutoff))==0 and len(threshold(all, cutoff))==0:
#                    if cutoff == 0.8:
#                        print "\tno significant mcm scores..."
#                    continue
#                print "\tProbality>=",cutoff, "yields", len(threshold(all,cutoff)),"mcm scores from",len(unique_ids(threshold(all,cutoff))),"unique sequences"
#                for scores,name in [(one,"best1"),(all,"all5")]:
#                    scores = threshold(scores,cutoff)
#                    good = correct(scores)
#                    accuracy =  prob(len(good),len(scores))
#                    #print scores
#                    #print good
#                    print "\t\t %s accuracy:"%name,accuracy
#        
#        # Now threshold them all
#        for cluster_id in sorted(unique_clusters(all_scores)):
#            cluster_scores = subset_clusters(all_scores, cluster_id)
#            
#            # Now look at the best scores
#            one = top_scores(cluster_scores, 1)
#            all = top_scores(cluster_scores)
#            sccs = [s for p,m,s,i,c in all][0]
#            
#            print "Cluster_id",cluster_id,"blast matches to",sccs,"with",len(unique_ids(cluster_scores)),"sequences"
#            for cutoff in cutoffs([8,9]):
#                if len(threshold(one, cutoff))==0 and len(threshold(all, cutoff))==0:
#                    if cutoff == 0.8:
#                        print "\tno significant mcm scores..."
#                    continue
#                print "\tProbality>=",cutoff, "yields", len(threshold(all,cutoff)),"mcm scores from",len(unique_ids(threshold(all,cutoff))),"unique sequences"
#                for scores,name in [(one,"best1"),(all,"all5")]:
#                    scores = threshold(scores,cutoff)
#                    good = correct(scores)
#                    accuracy =  div(len(good),len(scores))
#                    #print scores
#                    #print good
#                    print "\t\t %s accuracy:"%name,accuracy
 
        print "\n"
