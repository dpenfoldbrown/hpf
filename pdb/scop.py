# Defines a set of classes and functions for SCOP classifications
import string

class SCOPPart(object):
    """Defines a part of a scop domain on a single chain."""
        
    def __init__(self, chain, start=None, stop=None, chain_len=None):
        """
        Defines a portion of a scop definition on a single chain.  If the start
        or stop is undefined, the part is assumed to span the length of the
        chain.
        @param start: Location on the chain that starts the part.
        @param stop: Location on the chain that ends the part.
        @param chain_len: Length of the chain.  Used when start and stop are
            undefined.
        """
        self.chain = chain
        self.start = start
        self.stop = stop
        self.chain_len = chain_len
        #assert (self.start and self.stop) or self.chain_len
        if self.start or self.stop:
            assert self.start<self.stop
        
        #if self.stop == None or self.start == None:
        #    assert(self.stop == None and self.start == None)
        
    def length(self):
        """@return: The number of residues that comprise this part."""
        if not self.start or not self.start:
            assert(self.chain_len != None)
            return abs(self.chain_len)
        return abs(self.stop-self.start+1)

    def overlap(self, chain, start=None, stop=None):
        """
        See how much a start/stop alignment on a chain overlaps this part.
        If start or stop are none, your alignment is assumed to cover the entire
        chain.
        @param chain: The chain you aligned to.
        @param start: The start of your alignment to said pdb/chain.
        @param stop: The stop of your alignment to said pdb/chain.
        @return: The number of residues that overlap this scop classification.
        """
        if chain != self.chain:
            #print chain, self.chain
            return 0
        else:
            # If there is no start and stop, the overlap is the range queried.
            if not self.start or not self.start:
                return abs(start - stop)
            else:
                assert start < stop
                r=set(xrange(start,stop)).intersection(set(xrange(self.start,self.stop)))
                return len(r)
    
    def __repr__(self):
        if self.start or self.stop:
            end = str(self.start)+"-"+str(self.stop)
        else:
            end = ""
        return self.chain+":"+end
    
    def __eq__(self, other):
        return str(self) == str(other)

class SCOPDomain(object):
    """Defines a SCOP domain classification that may comprise of multiple parts
    across multiple chains."""
    
    def __init__(self, pdbid, sccs, parts=None, part_text=None,):
        """Create a SCOP domain object from a set of parts or text."""
        self.pdbid = pdbid
        self.sccs = sccs
        if parts:
            self.parts = parts
        else:
            try:
                self.parts = SCOPDomain._parse_text(part_text)
            except Exception, e:
                raise Exception(e,"Cannot parse text",part_text)
            
    def length(self):
        """Return the sum of the lengths of all parts"""
        return sum([part.length() for part in self.parts]) 
    
    @staticmethod
    def _parse_text(part_text):
        """Parses the part_text from ddb table."""
        # Examples: some parts are split up, chains in pdb are stupidly out of wack.
        # M:,N:
        # A:118-211,A:372-450
        # B:385-422,C:431-487
        p = []
        for part in part_text.split(","):
                parts = part.strip().split(":")
                chain = parts[0]
                # Now we work through the range part
                if len(parts[1]) > 1:
                    range_str = parts[1]
                    if range_str.startswith("-"):
                        range_str = range_str[1:]
                    try:
                        # Parse out junk
                        keep = string.digits+"-"
                        range_str = "".join(c for c in range_str if c in keep)
                        range_str = range_str.replace("S","")
                        start = int(range_str.split("-")[0])
                        stop = int(range_str.split("-")[1])
                    except:
                        print "error",part_text
                        raise
                # There's no range limit, use the whole chain
                else:
                    start = None
                    stop = None
                p.append(SCOPPart(chain, start,stop))
        return p
                        
    def overlap(self, pdbid, chain, start, stop):
        """Returns the number of residues overlapping this scop definition.
        @param pdbid: The pdb you aligned to.
        @param chain: The chain you aligned to.
        @param start: The start of your alignment to said pdb/chain.
        @param stop: The stop of your alignment to said pdb/chain.
        @return: The number of residues that overlap this scop classification.
        """
        assert(self.pdbid==pdbid)
        #print self.parts
        part_overlap = [part.overlap(chain,start,stop) for part in self.parts]
        #print part_overlap
        return sum(part_overlap)
    
    def __repr__(self):
        return self.pdbid+":"+self.sccs+"["+",".join(map(str, self.parts))+"]"
    
    def __eq__(self, other):
        return str(self) == str(other)
    
    
if __name__ == "__main__":
    print "Testing"
    domains = [SCOPDomain("pdb","sccs",part_text=part) for part in ["B:,A:","A:118-211,A:372-450","B:385-422,C:431-487"]]
    expected = {("A",0,400): [400-0, 211-118+400-372, 0],
               ("B",380,400): [400-380, 0, 400-385],
               ("D",385,422): [0,0,0]}
    for query in expected:
        chain,start,stop = query
        for i in xrange(len(expected[query])):
            try:
                result = domains[i].overlap("pdb",chain,start,stop)
                assert(expected[query][i] == result)
            except:
                print "query",query
                print "domain",domains[i]
                print "\tfound",result,"expected",expected[query][i]
                raise
    print "Finished"

