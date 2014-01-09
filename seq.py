'''
Created on Feb 24, 2010

@author: patrick
'''
from Bio.SeqIO import write
from tempfile import NamedTemporaryFile

class TemporaryRecordFile:
    """
    Write records to a temporary file that can be used in a "with" statement.
    """
        
    def __init__(self, records, format="fasta", *args, **kwargs):
        self.records = records
        self.format = format
        self.args = args
        self.kwargs = kwargs
    
    def __enter__(self):
        self.temp_file = NamedTemporaryFile(**self.kwargs)
        with open(self.temp_file.name,"w") as handle:
            write(self.records, handle, self.format)
        return self.temp_file.name
    
    def __exit__(self, type, value, traceback):
        return
        #self.temp_file.close()
