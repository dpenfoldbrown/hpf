#!/usr/bin/env python

# Example on how to query the uniprot ID translate web service

import urllib, urllib2

# For a full list of parameters, including all DB abbreviations:
# http://www.uniprot.org/faq/28
url = 'http://www.uniprot.org/mapping/'
params = {
    'from':'ACC+ID',
    'to':'HGNC_ID',
    'format':'tab',
    'query':''
}

data = urllib.urlencode(params)
request = urllib2.Request(url, data)
response = urllib2.urlopen(request)
page = response.read()
#line = response.readline()
#while line != '':
#    # Process line
#    line = response.readline()

print "Complete"

