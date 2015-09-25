#!/usr/bin/env python

# Unit tests for Psipred run functionality (used for MCM)
# dpb 7/26/2012

import unittest
from hpf.hddb.db import Session, Sequence, Psipred as PsipredORM

BASEPATH = "/Users/dpb/bonneau-dev/hpf/trunk/src/hpf/mcm/test/test_data/"

class PsipredGV4(unittest.TestCase):
    """Test a simple HPFPsipred run against the psipred prediction for that sequence already
    in the HPF database
    """
    def setUp(self, ):
        from hpf.pdb.psipred import HPFPsipredWrap
        session = Session()
        self.sequence_key  = 8560575
        self.ginzu_version = 4
        self.sequence = session.query(Sequence).get(self.sequence_key).sequence
        self.reference_pred = session.query(PsipredORM).filter_by(sequence_key=self.sequence_key, ginzu_version=self.ginzu_version).first().prediction

        self.db = BASEPATH + "small_db"

        # Run HPFPsipredWrap to get the new SS string
        self.created_pred = HPFPsipredWrap(sequence_key=self.sequence_key,
                                           nr_db=self.db,
                                           ginzu_version=self.ginzu_version,
                                           autorun=True,
                                           dbstore=False
                                          ).get_prediction_string()
    def tearDown(self, ):
        """Remove the files Psipred creates (.aux, .chk, .fasta, .mn, .mtx, .pn, .psipred, .psipred.1, .psipred.2, .sn)"""
        pass

    def testPrediction(self, ):
        """Check that the predicted secondary structure string is the same as the reference prediction (in the DB)"""
        print "Reference prediction: {0}".format(self.reference_pred)
        print "Created   prediction: {0}".format(self.created_pred)
        self.assertTrue(self.reference_pred == self.created_pred)

    def testPercentages(self, ):
        """Check the calculation of percent_alpha and _beta for the generated SS prediction string"""
        from hpf.pdb.psipred import percent_alpha_beta
        ref_length = 94
        ref_alpha = 0.074468085106382975
        ref_beta  = 0.48936170212765956

        length, alpha, beta = percent_alpha_beta(self.sequence, self.created_pred)

        print "[Field: Reference, Test]"
        print "Length: {0}, {1}".format(ref_length, length)
        print "%Alpha: {0}, {1}".format(ref_alpha, alpha)
        print "%Beta : {0}, {1}".format(ref_beta, beta)

        self.assertTrue(ref_length == length)
        self.assertTrue(ref_alpha == alpha)
        self.assertTrue(ref_beta == beta)

if __name__ == "__main__":
    unittest.main()

