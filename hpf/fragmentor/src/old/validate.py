import sys, os, re

def searchForPathExp(exp, path):
    dirList = os.listdir(path)
    regex = re.compile(exp)
    list = []
    for f in dirList :
        m = regex.match(f)
        if m != None:
            list.append(f)
    return list

def exists(path):
    return os.path.exists(path)

def validate(path):
    name = path.split("/")[-1]
    number = name[-3:]
    try:
        integer = int(number)
    except:
        return False, None, None
    
    aa = 'aa%s*' % name
    #print "Searching %s \n\tName: %s\n\tNumber:%s\n\tregex:%s" % (path, name, number, aa)
    list = searchForPathExp(aa, path)
    if len(list) == 2 :
        return True, name, integer
    else :
        return False, None, None
    

#aalz49503_05.075_v1_3
#aalz49509_05.075_v1_3

def main(argv):
    path=argv[0]
    
    dirList=os.listdir(path)
    finished = []
    for fname in dirList:
        #print fname
        valid, name, number = validate("%s/%s" % (path, fname))
        if valid :
            #print "%s is good" % path
            finished.append(number)
            
    completed = len(finished)
    finished.sort()
    max = finished[-1]
    min = finished[0]
    print "Completed: %d in range %d - %d" % (completed, min, max)
    span = range(min, max+1)
    unfinished = set(span)-set(finished)
    if len(unfinished) > 0 :
        l = []
        for i in unfinished :
            l.append(i)
        l.sort()
        print l
    
if __name__ == "__main__":
    main(sys.argv[1:])
