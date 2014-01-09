#!/usr/bin/env

# A script to create a matrix of astral-v-astral scores, in accordance
# with an astral index as created by create_astral_index.py. Stores a
# matrix (list of lists) in a pickle. Will enter mirror entries

import re
import cPickle as pickle

astral_list_file  = '/home/dpb/superfunc/astral/astral_mammoth.db.out'
astral_index_pkl  = '/home/dpb/superfunc/astral/astral_index.pkl'
astral_matrix_pkl = '/home/dpb/superfunc/astral/astral_matrix.pkl'
astral_matrix_columns_pkl = 'astral_matrix_columns.pkl'

# Open astral index pickle and load astral index and set matrix size
handle = open(astral_index_pkl)
astral_index = pickle.load(handle)
matrix_size = len(astral_index.keys())
handle.close()
    
def init_matrix():
    # Create and initialize astral matrix (fill with None objs)
    print "Initializing astral matrix"
    astral_matrix = []
    for i in range(matrix_size):
        column = []
        for j in range(matrix_size):
            column.append(None)
        astral_matrix.append(column)

    # Pickle empty astral matrix
    handle = open(astral_matrix_pkl, 'w')
    pickle.dump(astral_matrix, handle)
    handle.close()
    print "Matrix init complete"

def populate_matrix():
    print "Populating astral matrix"
    # Load empty matrix from pickle
    handle = open(astral_matrix_pkl)
    astral_matrix = pickle.load(handle)    
    handle.close()

    # Open astral list file
    handle = open(astral_list_file)
    
    # Add the score for each astral comparison to the matrix
    count = 0
    astral_list_pattern = r"(?P<src_astral>[a-zA-Z0-9_]+)\s(?P<cmp_astral>[a-zA-Z0-9_]+)\s(?P<score>[0-9.]+)"
    for line in handle:
        count += 1
        print "{0}\r".format(count), 
        found = re.match(astral_list_pattern, line)
        if not found:
            raise Exception("Line '{0}' doesn't match format".format(line))
        src_astral = found.group('src_astral')
        cmp_astral = found.group('cmp_astral')
        score = found.group('score')
    
        astral_matrix[astral_index[src_astral]][astral_index[cmp_astral]] = score
        #astral_matrix[astral_index[cmp_astral]][astral_index[src_astral]] = score
    
    handle.close()
    
    # Pickle astral matrix
    print "Pushing matrix columns to pickle file"
    handle = open(astral_matrix_columns_pkl, 'wb')
    for i in range(matrix_size):
        pickle.dump(astral_matrix[i], handle, pickle.HIGHEST_PROTOCOL)
        # Comment flush out to speed up
        #handle.flush()
        print "{0}\r".format(i), 
    handle.close()
    print "Matrix populate complete"

# Run init, let exit. Then change and run populate
#init_matrix()
populate_matrix()

