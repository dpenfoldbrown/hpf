import os
import subprocess
import numpy
from hpf.utilities import QueueParser
from tempfile import NamedTemporaryFile
from hpf.runtime import runtime

def interleave(file):
    from StringIO import StringIO
    modified = False
    with open(file) as handle:
        runtime().debug("interleave")
        lines = handle.readlines()
        if not lines[0].strip().split()[-1] == 'I':
            lines[0] = lines[0].replace("\n"," I\n")
            modified = True
        if modified:
            for line in lines[1:]:
                line = line.replace(" \n","\n")
                line = line.replace("X","-")
                line = line.replace("Z","-")
                
    if modified:
        with open(file,"w") as handle:
            map(handle.write,lines)
    

def clustal_to_phylip(src,dest=None):
    import os
    from Bio.AlignIO.ClustalIO import ClustalIterator
    from Bio.AlignIO.PhylipIO import PhylipWriter
    if dest==None:
        dest = src
    head,tail = os.path.split(dest)
    _dest = os.path.join(head,"."+tail)
    with open(src) as i:
        alignments = list(ClustalIterator(i))
    with open(_dest,"w") as f:
        writer = PhylipWriter(f)
        writer.write_file(alignments)
    os.rename(_dest,dest)
    return dest

class DnadistOptions(object):
    
    def __init__(self, input, output=None):
        self._itemp=None
        self.output = output
        if isinstance(input, str):
            assert(os.path.exists(input))
            self.input = input
        elif isinstance(input, list):
            from tempfile import NamedTemporaryFile
            from Bio.AlignIO.PhylipIO import PhylipWriter
            self._itemp = NamedTemporaryFile()
            self.input = self._itemp.name
            writer = PhylipWriter(self._itemp)
            writer.write_alignment(input)
            self._itemp.flush()
        else:
            raise Exception("Unknown input type",input)
        
    def stdin(self):
        from StringIO import StringIO
        input = StringIO()
        print >>input, "%s" % self.input
        print >>input, "P"
        print >>input, "P"
        #print >>input, "2"
        print >>input, "Y"
        return input.getvalue()
        
    def __del__(self):
        for temp in [self._itemp]:
            if temp:
                temp.close()

class Dnadist(object):
    
    def __init__(self, options):
        self.options = options
        
    def run(self):
        temp = NamedTemporaryFile()
        temp.write(self.options.stdin())
        temp.flush()
        cmd = "dnadist <%s >/dev/null" % (temp.name)
        print cmd
        subprocess.check_call(cmd, shell=True)
        temp.close()
        if self.options.output:
            os.rename("outfile",self.options.output)
        self.options.__del__()
        with open(self.options.output if self.options.output else "outfile") as array:
            return DistParser(array).parse().dict()

class ProtdistOptions(object):
    
    def __init__(self, input, output=None):
        self._itemp=None
        self.output = output
        if isinstance(input, str):
            assert(os.path.exists(input))
            self.input = input
        elif isinstance(input, list):
            from tempfile import NamedTemporaryFile
            from Bio.AlignIO.PhylipIO import PhylipWriter
            self._itemp = NamedTemporaryFile()
            self.input = self._itemp.name
            writer = PhylipWriter(self._itemp)
            writer.write_alignment(input)
            self._itemp.flush()
        else:
            raise Exception("Unknown input type",input)
        
    def stdin(self):
        from StringIO import StringIO
        input = StringIO()
        print >>input, "%s" % self.input
        print >>input, "P"
        print >>input, "P"
        #print >>input, "2"
        print >>input, "Y"
        return input.getvalue()
        
    def __del__(self):
        for temp in [self._itemp]:
            if temp:
                temp.close()

class Protdist(object):
    
    def __init__(self, options):
        self.options = options
        
    def run(self):
        temp = NamedTemporaryFile()
        temp.write(self.options.stdin())
        temp.flush()
        cmd = "protdist <%s >/dev/null" % (temp.name)
        print cmd
        subprocess.check_call(cmd, shell=True)
        temp.close()
        if self.options.output:
            os.rename("outfile",self.options.output)
        self.options.__del__()
        with open(self.options.output if self.options.output else "outfile") as array:
            return DistParser(array).parse().dict()
        
class DistParser(QueueParser):
    """Parses Protdist outfiles into a float array with a dictionary for matching names to rows/columns"""

    def __init__(self,handle):
        QueueParser.__init__(self,handle)
        self.n = int(self._get().strip())
        self.array = numpy.zeros((self.n,self.n))
        self.names = {}
        self.i = 0

    def next(self):
        if self._finished():
            raise StopIteration()
        name = self._peek().split()[0]
        self.names[self.i] = name
        j = 0
        while(j<self.n):
            line = self._get()
            for part in line.strip().split():
                try:
                    value = float(part)
                    self.array[self.i,j] = value
                    j+=1
                except ValueError:
                    assert(part==name)
        self.i += 1
        return (self.names[self.i-1],self.array[self.i-1])
        
    def parse(self):
        for res in self:
            pass
        return self
        
    def dict(self):
        d = {}
        for i in self.names:
            name1 = self.names[i]
            for j in self.names:
                name2 = self.names[j]
                d[(name1,name2)] = self.array[i,j]
        return d
        

class ProtparsOptions(object):
    
    def __init__(self, input, output=None):
        self._itemp=None
        self.output = output
        if isinstance(input, str):
            assert(os.path.exists(input))
            self.input = input
        elif isinstance(input, list):
            from tempfile import NamedTemporaryFile
            from Bio.AlignIO.PhylipIO import PhylipWriter
            self._itemp = NamedTemporaryFile()
            self.input = self._itemp.name
            writer = PhylipWriter(self._itemp)
            writer.write_alignment(input)
            self._itemp.flush()
        else:
            raise Exception("Unknown input type",input)
        
    def stdin(self):
        from StringIO import StringIO
        input = StringIO()
        print >>input, "%s" % self.input
        print >>input, "Y"
        return input.getvalue()
        
    def __del__(self):
        for temp in [self._itemp]:
            if temp:
                temp.close()

class Protpars(object):
    
    def __init__(self, options):
        self.options = options
        
    def run(self):
        temp = NamedTemporaryFile()
        temp.write(self.options.stdin())
        temp.flush()
        cmd = "protpars <%s >/dev/null" % (temp.name)
        print cmd
        subprocess.check_call(cmd, shell=True)
        temp.close()
        if self.options.output:
            os.rename("outtree",self.options.output)
        os.remove("outfile")
        self.options.__del__()

class ConsenseOptions(object):
    
    def __init__(self, input, output=None):
        self._itemp=None
        self.output = output
        if isinstance(input, str):
            assert(os.path.exists(input))
            self.input = input
        elif isinstance(input, list):
            from tempfile import NamedTemporaryFile
            from Bio.AlignIO.PhylipIO import PhylipWriter
            self._itemp = NamedTemporaryFile()
            self.input = self._itemp.name
            writer = PhylipWriter(self._itemp)
            writer.write_alignment(input)
            self._itemp.flush()
        else:
            raise Exception("Unknown input type",input)
        
    def stdin(self):
        from StringIO import StringIO
        input = StringIO()
        print >>input, "%s" % self.input
        print >>input, "Y"
        return input.getvalue()
        
    def __del__(self):
        for temp in [self._itemp]:
            if temp:
                temp.close()

class Consense(object):
    
    def __init__(self, options):
        self.options = options
        
    def run(self):
        temp = NamedTemporaryFile()
        temp.write(self.options.stdin())
        temp.flush()
        cmd = "consense <%s >/dev/null" % (temp.name)
        print cmd
        subprocess.check_call(cmd, shell=True)
        temp.close()
        if self.options.output:
            os.rename("outtree",self.options.output)
        os.remove("outfile")
        self.options.__del__()
