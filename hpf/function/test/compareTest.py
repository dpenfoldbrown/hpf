'''
Created on Mar 27, 2011

@author: slip
'''

import MySQLdb
import random

LLR_CUTOFF = 3

def main(*args):
    connection = _connection(db="hpf")
    connection2 = _connection(db="functionprediction")
    connection3 = _connection(db="hddb_IEA")
    global functionTermList
    
    global useTest
    useTest = True
    if useTest:
        testDict = _loadTestSet("/home/slip/Documents/MouseFunc_Key/go_testSet.txt")
    else:
        testDict = _loadTestSet("/home/slip/Documents/MouseFunc_Key/go_novelAnnotations.txt")
    
    print "Number of Genes in test file: " + str(len(testDict))
    
    correctCount = 0
    seqReturned = 0
    mfReturned = 0
    noFuncMatch = 0
    predictedButNotInTest = 0
    testAnnotations = 0
    numberOfPredsInDict = 0
    
    numberOfSeqWithHDDBTerms = 0
    numberOfHDDBMfs = 0;
    numberOfCorrectHddbMfs = 0;
    
    keyList = testDict.keys()
    switchList = list(keyList)
    random.shuffle(keyList)
    shuffDict = dict(zip(keyList,switchList))
    
    #print _fetchMFTermsFromSeqID("(133574L,)",connection3)
    
    

    for geneID in testDict.keys():
    #for geneID in 1:
        seqs= _fetchSeqFromGeneID(geneID,connection) # all sequence IDs attached to this gene ID
        foundPreds = False
        for seq in seqs:
            seqReturned = seqReturned + 1
            mfs = _fetchPredFromSeqID(str(seq),connection2,LLR_CUTOFF) # all func preds from this sequence
            realMfs = _fetchMFTermsFromSeqID(str(seq),connection3)
            
            #if len(mfs) > 0:
                #print "Gene: " + str(geneID) + "seq: " + str(seq)
                #print "mfs: " + str(mfs)
                #print "realMfs: " + str(realMfs)
            
            if len(realMfs) > 0:
                numberOfSeqWithHDDBTerms = numberOfSeqWithHDDBTerms + 1
                numberOfHDDBMfs += len(realMfs)
            #else:
                #print seq

            for mf in mfs:
                if not(foundPreds):
                    foundPreds = True
                    testAnnotations = testAnnotations + len(testDict[geneID]) # number of annotations on genes where we found a func pred
                mfReturned = mfReturned + 1
                mf = str(mf)[2:-3]
                if functionTermList.__contains__(mf):
                    numberOfPredsInDict = numberOfPredsInDict + 1
                if testDict[geneID].__contains__(mf): # if true we have a match
                #if testDict[shuffDict[geneID]].__contains__(mf):
                    correctCount = correctCount + 1
                    if mf == "GO:0005554": # This would be a match to unknown bio process
                        noFuncMatch = noFuncMatch + 1
                else:
                    predictedButNotInTest = predictedButNotInTest + 1; # our prediction is not in the test set annotation
                if realMfs.__contains__(mf):
                    numberOfCorrectHddbMfs +=1
                    
    
    
    print "Number of sequences returned from test gene ids: " + str(seqReturned)
    print "Number of mf predictions returned from sequence keys: " + str(mfReturned)
    print "Number of func preds with corresponding column: " + str(numberOfPredsInDict)
    #print "Number of matches that were unknown function: " + str(noFuncMatch)
    print "Number of correct predictions: " + str(correctCount)
    print "Number of func preds not in the test set annotations: " + str(predictedButNotInTest)
    print "Number of total test set annotations for genes with returned predictions: " + str(testAnnotations)
    print "Number of sequences with Hddb entries: " + str(numberOfSeqWithHDDBTerms)
    print "Number of Hddb mf term matches: " + str(numberOfCorrectHddbMfs)
    print "Number of Hddb mf terms: " + str(numberOfHDDBMfs)

    
    
# function to load the test data from a file 
def _loadTestSet(testSetLoc):
    print testSetLoc    
    f = open(testSetLoc,'r')
    global functionTermList
    functionTermList = f.readline()
    functionTermList = functionTermList.split()
    for i in range(1,len(functionTermList)):
        functionTermList[i]= functionTermList[i].split('.')[0][1:]
    
    testGeneDict = dict()
    global useTest
    for line in f:
        line = line.split()
        if not(useTest):
            line[0] = line[0].split("\"")[1] # this line is necessary in novelAnnotations because ids have quotes
        testGeneDict[line[0]] = [];
        for entryNum in range(1,len(line)):
            if line[entryNum]=="1":
                testGeneDict[line[0]].append(functionTermList[entryNum])
    
    f.close()
    
    return testGeneDict
             

    
# function to get seq key from gene id    
def _fetchSeqFromGeneID(geneID, connection):

    cursor = connection.cursor()
    query = "select distinct(sequence_key) from mouseIDmap where gene_id=\"" + geneID + "\""
    #print query
    cursor.execute(query)
    return cursor.fetchall()

# function to get mf predictions from seq key
def _fetchPredFromSeqID(seqID,connection,llrCutoff):
    cursor = connection.cursor()
    query = "select mf_acc from bayes_golite_062009_3 where parent_sequence_key=" + seqID[1:-3] + " and pls_llr>" +str(llrCutoff) + " and type=\"psi\""
    #print query 
    cursor.execute(query)
    return cursor.fetchall()

# function to get actual mf annotations from out db for a sequence key
def _fetchMFTermsFromSeqID(seqID,connection):
    cursor = connection.cursor()
    query = "select acc from hddb_IEA_goLite_082010_medconf where sequence_key=" + seqID[1:-3] + " and term_type=\"molecular_function\""
    #print query
    cursor.execute(query)
    realTerms = list()
    for term in cursor.fetchall():
        #print term
        realTerms.append(str(term)[2:-3])
    return realTerms

    
def _connection(db=None):
    if db:
        _db = MySQLdb.connect(db=db,read_default_file="~/.my.cnf")
    else:
        _db = MySQLdb.connect(read_default_file="~/.my.cnf")
        #_db = MySQLdb.connect(host="127.0.0.1", user="pwinters", passwd="bonneaulab", port=13307,db="hpf")
    return _db
   
def _print(*args):  
    print " ".join([str(a) for a in args])
 
    
if __name__=="__main__":
    main()

