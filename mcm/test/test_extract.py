#!/usr/bin/env python

# Unit tests for Rosetta extract functionality
# dpb 7/25/2012

import os
import unittest
from hpf.mcm.extract import extract

BASEPATH = "/home/dpb/mcm/mcm/test/test_data/"

class Rosetta2Extractor(unittest.TestCase):
    """Will test the accuracy of the Rosetta++ 2.? extractor function from hpf.mcm.extract
    Note: Must have rosetta++ built and executable, and have the rosetta++ database directory.
    """
    def setUp(self, ):
        self.reference_pdb = BASEPATH + 'S_0003_9566.pdb.reference'

        self.silent_file   = BASEPATH + 'na977.top100'
        self.tag        = 'S_0003_9566'
        self.log_file   = 'extract.log'
        self.paths_file = '/home/dpb/mcm/env/paths.txt'

        # Run extract on test files
        self.created_pdb = extract(source_file=self.silent_file, tag=self.tag, log_file=self.log_file, paths=self.paths_file, debug=True)

    def tearDown(self, ):
        """Clear the log and pdb files created to clean for the next test"""
        import os
        print "Cleaning up after test.."
        os.remove(self.log_file)
        os.remove(self.created_pdb)
    
    def testFileExistence(self, ):
        """Check that log and PDB file were created"""
        import os
        self.assertTrue(os.path.isfile(self.log_file))
        self.assertTrue(os.path.isfile(self.created_pdb))

    def testFileIdentity(self, ):
        """Check that both files are identical to reference files.
        Note that this checks to see if the files are EXACTLY the same.
        There is a better way to do this.
        """
        ref_handle = open(self.reference_pdb)
        new_handle = open(self.created_pdb)
        
        ref_lines = ref_handle.read().split('\n')
        new_lines = new_handle.read().split('\n')

        self.assertTrue(len(ref_lines) == len(new_lines))

        match = True
        for i in range(len(ref_lines)):
            if ref_lines[i] != new_lines[i]:
                match = False
                break
        self.assertTrue(match)


if __name__ == "__main__":
    unittest.main()
