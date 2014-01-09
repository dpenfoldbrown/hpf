'''
Created on Apr 28, 2010

@author: Patrick
'''
from tempfile import NamedTemporaryFile

class ControlFile(object):

    def __init__(self,file,**kwargs):
        self.file = file
        self.temp = False if self.file else True
        self.kwargs = kwargs
    
    def _write(self,handle):
        assert False, "Abstract method"
    
    def __enter__(self):
        if self.temp:
            self.temp_file = NamedTemporaryFile(**self.kwargs)
            self.file = self.temp_file.name
        with open(self.file,"w") as handle:
            self._write(handle)
        return self.file
    
    def __exit__(self, type, value, traceback):
        if self.temp:
            self.temp_file.close()