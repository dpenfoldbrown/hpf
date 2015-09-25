#!/usr/bin/env python
import unittest
from hpf.function.scripts import ypdr_format

class test_ypdr_format(unittest.TestCase):
    
    def test_format(self):        
        euk_fr_spec = {"type":"fr",
               "experiment_key":804,
               "base_llr":-15,
               "pls_llr":12,
               "parent_sequence_key":0,
               "domain_sequence_key":1,
               "name":"molecular function",
               "acc":"GO:000000",
               "timestamp":"NOW"}
        result = ypdr_format.format(euk_fr_spec)
        self.failUnlessEqual(result["llr"], 0, "Should be adjusted to 0: "+str(result["llr"]))
        
        prok_mcm_gen = {"type":"mcm",
               "experiment_key":809,
               "base_llr":-1,
               "pls_llr":5.25,
               "parent_sequence_key":0,
               "domain_sequence_key":1,
               "name":"molecular function",
               "acc":"GO:000000",
               "timestamp":"NOW"}
        result = ypdr_format.format(prok_mcm_gen)
        self.failUnlessEqual(result["llr"], 0, "Should be adjusted to 0: "+str(result["llr"]))
        
if __name__ == "__main__":
    unittest.main()