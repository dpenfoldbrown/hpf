'''
Created on Feb 24, 2010

@author: patrick
'''

from Bio.Align.Generic import Alignment as GenericAlignment
from Bio import AlignIO
from tempfile import NamedTemporaryFile

class TemporaryAlignmentFile:
    """
    Write alignments to a temporary file that can be used in a "with" statement.
    """
    
    def __init__(self, alignments, format="fasta", *args, **kwargs):
        self.alignments = alignments
        self.format = format
        self.args = args
        self.kwargs = kwargs
    
    def __enter__(self):
        self.temp_file = NamedTemporaryFile(**self.kwargs)
        with open(self.temp_file.name,"w") as handle:
            format = "phylip" if self.format=="phylipi" else self.format
            AlignIO.write(self.alignments, handle, format)
        if self.format == "phylipi":
            from hpf.phylip import interleave
            interleave(self.temp_file.name)
        return self.temp_file.name
    
    def __exit__(self, type, value, traceback):
        return
        #self.temp_file.close()


class AlignmentFactory(object):
    
    def __init__(self, class_):
        self._class = class_
    
    def copy(self, alignment):
        _copy = self._class(alignment._alphabet)
        if hasattr(_copy,"_copy"):
            _copy._copy(alignment)
        else:
            _copy._records = alignment._records
        return _copy
    
    def read(self, handle, format="fasta", **kwargs):
        return self.copy(AlignIO.read(handle, format, **kwargs))
    
    def parse(self, handle, format="fasta", **kwargs):
        return map(self.copy,AlignIO.parse(handle, format, **kwargs))

class Alignment(GenericAlignment, object):
    def __init__(self, *args, **kwargs):
        GenericAlignment.__init__(self, *args, **kwargs)
        
    def _copy(self, alignment):
        self._records = alignment._records 
