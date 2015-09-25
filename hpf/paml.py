import os
from Queue import Queue
from tempfile import NamedTemporaryFile
from hpf.runtime import runtime

def codeml_label(tree,prefix="node"):
    """
    Names nodes 'taxon' according to the PAML CodeML strategy.  
    Each node.data.taxon=prefix+CodeMLID
    """
    # Now we need to figure out ancestor names and id's for CodeML
    # CodeML uses a different numbering system for nodes where it
    # numbers from the sequence file first, then incrementally over
    # nodes using a BFS
    queue = Queue()
    queue.put(tree.root)
    # Start naming at number of sequences+1 like codeML
    count = len(tree.get_terminals())+1
    while not queue.empty():
        id = queue.get()
        node = tree.node(id)
        if not node.data.taxon:
            node.data.taxon = prefix+str(count)
            count += 1
            for succ in node.get_succ():
                queue.put(succ)
    return tree

class CodeMLParser(object):
    
    def expression(self):
        from pyparsing import Suppress,Combine,Optional,oneOf,OneOrMore,Word,nums,Group,alphas,alphanums,Literal,SkipTo,empty,lineEnd
        cvtInt = lambda toks: int(toks[0])
        cvtReal = lambda toks: float(toks[0])
        cvtTuple = lambda toks : tuple(toks.asList())
        nameJoin = lambda toks : "".join([tok.replace("#","") for tok in toks[0]])
        #lambda toks: " ".join([str(t) for t in toks[0]])
        
        # define punctuation as suppressed literals
        lparen,rparen,lbrack,rbrack,lbrace,rbrace,colon = map(Suppress,"()[]{}:")

        integer = Combine(Optional(oneOf("+ -")) + Word(nums))\
            .setName("integer")\
            .setParseAction( cvtInt )
        real = Combine(Optional(oneOf("+ -")) + Word(nums) + "." +
                       Optional(Word(nums)) +
                       Optional(oneOf("e E")+Optional(oneOf("+ -")) +Word(nums))).setName("real").setParseAction( cvtReal )

        # TREE DEFINITION
        # ((seq2: 0.537243, seq1: 0.000004): 0.255741, seq3: 0.281503);
        tree_w_branches = (
            OneOrMore(Word("():,."+alphas+nums))+Literal(";")
            ).setParseAction(lambda tokens: " ".join(tokens[:-1])+";")

        # SITE PROBABILITIES
        # site Freq Data:
        # 1    1    AAA: A(0.978) A(1.000)
        site_prob = (
            integer.setResultsName("site",listAllMatches=True) + 
            integer.setResultsName("freq",listAllMatches=True) + 
            Word(alphas+"-").setResultsName("extant",listAllMatches=True) + colon + 
            Group(OneOrMore(Group(Word(alphas,exact=1)+lparen+real+rparen))).setResultsName("probability",listAllMatches=True) +
            lineEnd
            )

        # ANCESTRAL SEQUENCES
        # seq1       ACC
        # node #4    ACC
        # Optional # character with node # needs to be joined into a single name
        sequence =  (
                     Group(Word(alphanums)+
                           Optional(Combine(Literal("#")+Word(nums)))).setParseAction(nameJoin).setResultsName("name",listAllMatches=True)+
                           Word(alphas+"- ").setResultsName("sequence", listAllMatches=True)+lineEnd
                    )
                     
        
        return (SkipTo(Literal("Ancestral reconstruction by AAML."),include=True).suppress() +
                tree_w_branches.setResultsName("tree") +
                SkipTo(Literal("site")+Literal("Freq")+Literal("Data:"), include=True,).suppress()+
                Group(OneOrMore(site_prob)).setResultsName("sites")+
                SkipTo(Literal("List of extant and reconstructed sequences")+Word(nums)+Word(nums), include=True).suppress()+
                Group(OneOrMore(sequence)).setResultsName("sequences")+
                SkipTo(Literal("for a site."),include=True).suppress()+
                Group(OneOrMore(real)).setResultsName("probability")+
                empty
                )
        
        
        
    def parse(self, input):
        """
        @return: (alignment,probability,site_mat)
            site_mat contains probability of site for each ancestral sequence
        """
        from treecorr.tree import TreeAlignment
        from Bio.Alphabet import generic_protein
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        import numpy
        result = self.expression().parseFile(input)
        
        tree_str = result['tree']
        num_ancestral = len(result['probability']) 
        num_extant = len(result['sequences']['name'])-num_ancestral
        alignment = TreeAlignment(tree_str,alphabet=generic_protein)
        for i,name in enumerate(result['sequences']['name']):
            alignment.add_sequence(name,result['sequences']['sequence'][i].replace(" ",""))
            record = alignment._records[-1]
            record.annotations["probability"] = 1.0 if i<num_extant else result['probability'][i-num_extant]
            
            # Add in site probabilities to the sequence record. 
            if i>=num_extant:
                site_prob = numpy.zeros(len(record.seq),dtype=float)
                for site,r in enumerate(result['sites']['probability']):
                    node = i-num_extant
                    aa,prob = r[node]
                    site_prob[site] = prob
            else:
                site_prob = numpy.array([1.0 for i in xrange(len(record.seq))],dtype=float)
            
            if not hasattr(record, "letter_annotations"):
                record.letter_annotations = {}
            record.letter_annotations["probability"] = site_prob        
                    
        return codeml_label(alignment)
        #print num_ancestral, num_extant 
        #probability = [1.0 if i<num_extant else result['probability'][i-num_extant] for i in xrange(num_extant+num_ancestral)]
        #site_mat = numpy.zeros((num_ancestral,len(result['sites']['site'])))
        #for site in result['sites']['site']:
        #    for node,r in enumerate(result['sites']['probability'][site-1]):
        #        aa,prob = r
        #        print node,site,aa,prob
        #        site_mat[node,site-1] = prob
        #return (alignment, site_mat)
        
class ControlFile(object):
    
    def __init__(self,base_control_file,options=None,file=None,**kwargs):
        self.base_control_file = base_control_file
        self._parse(options)
        self._file = file
        self._kwargs = kwargs

    def _parse(self,options):
        self.options={}
        with open(self.base_control_file) as handle:
            for line in handle:
                line = line.strip()
                line = line[0:line.index("*")] if line.find("*") != -1 else line
                parts = line.split("=")
                if len(parts)==2:
                    name = parts[0].strip()
                    value = parts[1].strip()
                    self.options[name]=value
        if options:
            self.options.update(options)
    
    def __enter__(self):
        if not self._file:
            self._temp_file = NamedTemporaryFile(**self._kwargs)
            self._file = self._temp_file.name
        with open(self._file,"w") as handle:
            for key in self.options:
                print >>handle, "%s = %s" % (key,self.options[key])
        return self._file
    
    def __exit__(self, type, value, traceback):
        return
    
        
class CodeML(object):

    def __init__(self, ctl, dir=None):
        self.ctl = ctl
        if dir==None:
            from tempfile import mkdtemp
            self.dir = mkdtemp()
        else:
            self.dir = dir
    
    def run(self):
        import subprocess,os
        with self.ctl as ctl_file:
            cmd ="codeml "+ctl_file
            runtime().debug(cmd)
            subprocess.check_call(cmd,shell=True,cwd=self.dir,stdout=open(os.devnull,"w"))
        with open(os.path.join(self.dir,"rst")) as handle:
            return CodeMLParser().parse(handle)

class RSTParser(object):
    """Parses RST file for ancestral sequences."""
        
    def __init__(self, file_handle, prefix="node"):
        """
        @deprecated: Use RSTParser
        """
        self.file_handle = file_handle
        self.ids = None
        self.prefix = prefix
        __start = "List of extant and reconstructed sequences"
        line=""; start=False
        while not start:
            line = self.file_handle.readline().strip()
            if line.startswith(__start):
                #print "Good",line
                start = True
            else:
                pass
            #print "Bad",line
        if not start:
            raise RuntimeError("Cannot read RST file correctly")
        self.__start()
       
    def __iter__(self):
        """@return: self._seqs()"""
        return self._seqs()

    #def alignment(self):
    #    alignment = TreeAlignment()
    #    for name,seq in self:
    #        alignment.add_sequence(name,seq)
    #    return alignment
        

    def __start(self):
        """
        Expects a well defined RST file.
        """
        self.file_handle.readline()
        line = self.file_handle.readline()
        number, length = line.strip().split()
        self.number = int(number)
        self.length = int(length)
        self.file_handle.readline()
        #print number,length
        
    def _seqs(self):
        """Generator yields all (name,seq) pairs for an RST file"""
        number = self.number
        for i in xrange(self.number):
            line = self.file_handle.readline().strip()
            number -= 1
            yield(self._seq(line))
        if number != 0:
            raise RuntimeError("Couldn't find the number of expected sequences.")
        
    def _seq(self,line):
        """Returns (name,seq) pair for a line"""
        name = line[:15].strip()
        if name.startswith("node"):
            name = self.prefix+name[6:]
        seq = "".join(line[15:].strip().split())
        #print name,seq
        assert(len(seq) == self.length)
        return (name,seq)

    def alignment(self, alignment):
        """
        Add the sequences in this parser to an alignment.
        @return: The alignment
        """
        for name,seq in self:
            alignment.add_sequence(name,seq)
        return alignment
