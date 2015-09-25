'''
Created on May 12, 2010

@author: Patrick
'''
import unittest
from hpf.hddb.mcm.mammoth import MammothScore
from hpf.hddb.mcm import McmDB,McmFactory
from hpf.hddb.mcm.cluster import ClusterCenter

class TestMcmFactory(unittest.TestCase):
#                        id: 3277636
#              sequence_key: 8555583
#               outfile_key: 240882
#    outfile_mcm_result_key: NULL
#                      scop: 1.69
#             structure_key: 1115640
#          prediction_index: 2
#       n_decoys_in_outfile: 10000
#       cluster_center_rank: 0
#   experiment_percent_beta: 0
#                    target: mz899.result.15000.10000.out
#                      psi1: 96.4
#  prediction_contact_order: 4.1
#                      psi2: 42
#      experiment_astral_ac: d1jt6a2
#                    bratio: 148.41
#               convergence: 9.206445
#           prediction_file: decoy_1876.pdb
#                decoy_name: decoy_1876.pdb
#  prediction_percent_alpha: 0.132743362831858
#                    evalue: 0.007096
#                     class: 2
#          experiment_index: 2898
#   experiment_sequence_key: 14791
#experiment_sequence_length: 115
#  experiment_percent_alpha: 0.747826086956522
#                      ln_e: 4.948
#           experiment_sccs: a.121.1.1
#prediction_sequence_length: 112
#                    aratio: 0.177505659600741
#                    zscore: 4.831
#                     ratio: 0.973913043478261
#       cluster_center_size: 200
#                       nss: 36
#               probability: 0.191084454720235
#                     score: 560.2
#  experiment_contact_order: 4.5
#   prediction_percent_beta: 0.283185840707965
#           experiment_file: 18909.pdb
#                      nsup: 108
#      cluster_center_index: 1876
#                 timestamp: 0000-00-00 00:00:00
#                       tag: 0

# Experiment, mcmdb

#*************************** 1. row ***************************
    dbentry = {
                  'percent_beta': 0.0,
                     'astral_ac': 'd1jt6a2',
                  'sequence_key': 0,
                        'length': 115,
                 'percent_alpha': 0.747826086956522,
                          'sccs': 'a.121.1.1'
               }


# Prediction
    prediction = {
                  'sequence_key': 0,
                 'structure_key': 0,
           'n_decoys_in_outfile': 0,
           'cluster_center_rank': 0,
                        'target': '',
                   'convergence': 9.206445,
      'prediction_percent_alpha': 0.132743362831858,
    'prediction_sequence_length': 112,
           'cluster_center_size': 0,
       'prediction_percent_beta': 0.283185840707965,
          'cluster_center_index': 0,
                   'outfile_key': 0,
               'prediction_file': 'decoy_1876.pdb'

                   }

# Mammoth
    score = { 'prediction_index': 0,
                          'psi1': 0,
      'prediction_contact_order': 4.1,
                          'psi2': 0,
                    'decoy_name': 'decoy_1876.pdb',
                        'evalue': 0,
              'experiment_index': 0,
                          'ln_e': 0,
                        'zscore': 4.831,
                           'nss': 0,
                         'score': 0,
      'experiment_contact_order': 0,
                          'nsup': 0,
                    'experiment': '1.pdb',
                    'prediction': 'decoy_1876.pdb'
            }
# MCM
#                     class: 3
#                    aratio: 0.654329147388994
#                    bratio: 148.41
#                     ratio: 0.888888888888889
#               probability: 0.191084454720235
#           experiment_file: 987422.pdb


    def setUp(self):
        self.mcmdb = McmDB()
        self.mcmdb.add(1,**TestMcmFactory.dbentry)
        self.convergence = {}
        from hpf.hddb.db import McmData
        self.convergence[TestMcmFactory.score['decoy_name']] = ClusterCenter.create(McmData(**TestMcmFactory.prediction))

    def tearDown(self):
        pass


    def testProbability(self):
        score = MammothScore(**TestMcmFactory.score)
        factory = McmFactory(self.convergence, self.mcmdb)
        data = factory.create(score)
        assert data.__dict__["class"] == 2, data.__dict__["class"]
        # allowing for numerical differences between the different languages. 
        assert abs(data.probability-0.191084454720235)<1e-7, str(data.probability)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()