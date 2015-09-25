#!/usr/bin/env python

# Creates a dictionary of the form (astral_sid -> index int), where the index of
# the astral_sid is the position of that astral_sid in a string-sorted list of
# unique astral_sids.
# Dictionary is then stored in a pickle.
# For use with create_astral_matrix and ultimately cluster_struct_allvall.py
# to create flat files so astral all-v-all info can be read into memory

import re
import cPickle as pickle


astral_list_file = '/Users/dpb/Documents/superfunc/data/humanrna/astral_mammoth.db.out'
astral_index_pkl = 'astral_index.pkl'

# Open astral list file
handle = open(astral_list_file)

# Create set from astral IDs in list
astral_list_pattern = r"(?P<src_astral>[a-zA-Z0-9_]+)\s(?P<cmp_astral>[a-zA-Z0-9_]+)\s[0-9.]+"
unique_astrals = set()
for line in handle:
    found = re.match(astral_list_pattern, line)
    if not found:
        raise Exception("Line '{0}' doesn't match format".format(line))
    unique_astrals.add(found.group('src_astral'))
    unique_astrals.add(found.group('cmp_astral'))

handle.close()
print "Number of unique astral IDs: {0}".format(len(unique_astrals))

# Create a sorted list from the unique astrals
sorted_astrals = sorted(unique_astrals)

# Add the astrals and their indexes in sorted order to a dictionary
astral_index = {}
for i in range(len(sorted_astrals)):
    astral_index[sorted_astrals[i]] = i

print "Number of keys in astral index: {0}".format(len(astral_index.keys()))

# Output astral index to pickle
handle = open(astral_index_pkl, 'w')
pickle.dump(astral_index, handle)
handle.close()

print "Complete"


