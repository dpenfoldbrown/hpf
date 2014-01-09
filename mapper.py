class MapperInterface(object):
    """
    Base interface for defining index mappers.  These are used to do things such
    as map domain indices to alignment columns, pdb sequence records to atom records, etc.
    Can be used as a dictionary, inverted, or iterated.

    Define the "mapping" method to create a mapper and make sure that it always yields a
    forward mapping.  Inverting is handled by the __iter__ and __getitem__ methods.
    """

    def __init__(self, inverse=False):
        """
        @param inverse: 
        """
        self._forward = None
        self._reverse = None
        self._inverse = inverse

    def __iter__(self):
        """
        Iterate over the mapping tuples, inverting their direction if specified.
        Yields (a_index,b_index) tuples or (b_index,a_index) if inverted.
        No caching is performed.
        """
        for t in self.mapping():
            yield tuple(reversed(t)) if self._inverse else t

    def mapping(self):
        """
        Yields a (a_index,b_index) tuple iterable that maps 
        indices from A to indices in B.  This is never inverted, but defines the
        standard method for determining the mapping.
        """

    def inverse(self):
        """
        Returns a shallow copy of this mapper that will perform lookups or iterate
        in the reverse direction.  Sets self._inverse=True on the returned copy.
        """
        from copy import copy
        z = copy(self)
        z._inverse = True
        return z
        
    def _dict_(self):
        """
        Generate the cached dictionary mapping
        """
        if self._forward==None or self._reverse==None:
            self._forward = {}
            self._reverse = {}
            for a,b in self.mapping():
                self._forward[a]=b
                self._reverse[b]=a
        return self._reverse if self._inverse else self._forward

    def __getitem__(self, a):
        d = self._dict_()
        return d[a]

class MapperChain(MapperInterface):
    """
    Index mappings from one sequence to another are useful when mapping, for example,
    PDB sequence records to atom records or domain indices to protein indices.
    This class accepts some list of mappers and chains them together.
    """

    def __init__(self, *mappers, **kwargs):
        """
        @param mappers: A variable number of objects of the mapper interface that will
        be chained in order to create a mapping from mappers[0]-> mappers[1]->...->mappers[-1]
        @type mappers: list<MapperInterface>
        """
        super(MapperChain,self).__init__(**kwargs)
        self._mappers = mappers

    def _dict_(self):
        if self._forward==None or self._reverse==None:
            pred = None # Predescessor lookup
            for mapper in self._mappers:
                temp = {}
                # We get to z though a's predescessor
                for a,z in mapper:
                    #print a,z
                    if pred==None or pred.has_key(a):
                        # AA is a's predescessor
                        AA = pred[a] if pred else a
                        temp[z] = AA 
                    else:
                        #print "pred does not have key",a
                        pass
                #print "round",temp
                pred = temp

            # Store the inverse of the above as the feed-forward dictionary.
            #print "done,saving",pred
            self._forward = {}
            for end in pred:
                #print "Finished",pred[end],end
                self._forward[pred[end]]=end
            self._reverse = pred
        return self._reverse if self._inverse else self._forward

    def mapping(self):
        # Ensure the dictionary has been calculated
        self._dict_()
        # Use the forward mapper anyway
        d = self._forward
        for key in d:
            yield (key,d[key])   

class ProteinDomainMapper(MapperInterface):
    """
    Uses a domain's region object to map protein indices to
    domain indices.
    """
    def __init__(self, domain, **kwargs):
        """
        @type domain: hpf.hddb.db.Domain
        """
        super(ProteinDomainMapper,self).__init__(**kwargs)
        self.domain = domain

    def mapping(self):
        for domain,protein in enumerate(xrange(self.domain.region.start,self.domain.region.stop)):
            yield (protein,domain)

class ProteinAlignmentMapper(MapperInterface):
    """
    Maps protein indices to an alignment's columns.
    """

    def __init__(self, protein, alignment, **kwargs):
        """
        @type alignment: Bio.Align.Generic.Alignment
        @type protein: hpf.hddb.db.Protein
        @type domain: hpf.hddb.db.Domain
        """
        super(ProteinAlignmentMapper,self).__init__(**kwargs)
        self._protein = protein
        self._alignment = alignment
        self._alignment_index = None
        # Find the protein in the alignment
        # for i,seq in enumerate(alignment):
        #     print seq.id, self._protein.id

        for i,seq in enumerate(alignment):
            if seq.id==str(self._protein.id):
                self._alignment_index = i
        assert self._alignment_index != None


    def mapping(self):
        """
        Yields a (residue,column) tuple generator that maps 
        residues to columns in the alignment.
        """
        protein_record = self._alignment[self._alignment_index]
        protein_seq = str(protein_record.seq)
        aa_count=0
        for column,aa in enumerate(protein_seq):
            if aa!='-':
                aa_count=aa_count+1
                yield(aa_count,column)

class DomainAlignmentMapper(ProteinAlignmentMapper):
    """
    Short circuit for mapping a Domain to an alignment.  Uses the
    ProteinAlignmentMapper and only yields domain residues indexed at 0.
    """
    def __init__(self, protein, domain, alignment, **kwargs):
        super(DomainAlignmentMapper,self).__init__(protein,alignment, **kwargs)
        self._domain = domain

    def mapping(self):
        domain_region = self._domain.region
        for res,column in super(DomainAlignmentMapper,self).mapping():
            if res>=domain_region.start and res<=domain_region.stop:
                yield (res-domain_region.start,column)

    

if __name__ == "__main__":
    # Test some stuff
    class TestMapper(MapperInterface):
        def __init__(self, *tuples, **kwargs):
            super(TestMapper,self).__init__(**kwargs)
            self.tuples = tuples
        def mapping(self):
            return iter(self.tuples)
        
    m1 = TestMapper(('A','a'),('B','b'),('C','c'),('E','e'))
    m3 = TestMapper(('A','a'),('B','b'),('C','c'),('E','e'), inverse=True)
    m2 = TestMapper(('d','D'),('c','Z'),('a','A'))
    chain = MapperChain(m1,m2)
    chain._dict_()
    print "made chain"
    inv = chain.inverse()
    print "Forward",chain._forward
    print "Inverse",inv._forward
    assert chain['A']=='A'
    assert chain['C']=='Z'
    assert inv['Z']=='C'
    for a,b in chain:
        print (a,b)
    assert list(m1)!=list(m2)
    assert list(m1.inverse())==list(m3)
