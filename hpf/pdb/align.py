'''
Created on Jan 27, 2010

@author: Patrick
'''
from hpf.pdb import get_seq

class IdentityMatrix():
    pass

class StructureAlignment(object):
    
    def __init__(self, prediction, experiment, mapping):
        self.prediction = prediction
        self.experiment = experiment
        self.mapping = mapping
        self.pred_seq = get_seq(self.prediction)
        self.exp_seq = get_seq(self.experiment)
        
    def conservation(self, inspect_sites=None, matrix=None, gap= -11, neighbor=0):
        """
        Determine sequence conservation based on the structure alignment.
        @param inspect_sites: sites on the experiment to inspect, otherwise all. 
        @return: generator len(mapping), matrix values or True/False identical and prediction index.
        """
        pred_seq = self.pred_seq
        exp_seq = self.exp_seq
        
        def neighbor_scores(pred_site, exp_site):
            """Examine close neighbors on the experiment sequence for conservation"""
            #for pred in range(pred_site-neighbor,pred_site+neighbor+1):
            for exp in range(exp_site-neighbor,exp_site+neighbor+1):
                if exp>=0 and exp<len(exp_seq):#pred>=0 and pred<len(pred_seq) and 
                        yield (matrix(pred_seq[pred_site], exp_seq[exp]),pred_site,exp_site) if matrix else (pred_seq[pred_site] == exp_seq[exp], pred_site,exp_site)
        
        inspect_sites = inspect_sites if inspect_sites else xrange(len(exp_seq))
        #print "inspecting",len(inspect_sites)
        for pred_site, exp_site in self.mapping:
            if exp_site in inspect_sites:
                if pred_site == None or exp_site == None:
                    yield (gap,pred_site,exp_site) if matrix else (False, pred_site, exp_site)
                else:
                    yield max(neighbor_scores(pred_site, exp_site))
                    #yield (matrix(pred_seq[pred_site], exp_seq[exp_site]),pred_site) if matrix else (pred_seq[pred_site] == exp_seq[exp_site], pred_site)
                

                
