'''
Created on Mar 16, 2010

@author: patrick
'''
from tempfile import NamedTemporaryFile

class TemporaryTreeFile:
    """
    Write alignments to a temporary file that can be used in a "with" statement.
    """
    
    def __init__(self, tree, format="newick", *args, **kwargs):
        self.tree = tree
        self.format = format
        self.args = args
        self.kwargs = kwargs
    
    def __enter__(self):
        self.temp_file = NamedTemporaryFile(**self.kwargs)
        with open(self.temp_file.name,"w") as handle:
            if self.format == "newick":
                io = self.tree.to_string(plain_newick=True)+";"
            else:
                io = self.tree.to_string()
            print >>handle, io
        return self.temp_file.name
    
    def __exit__(self, type, value, traceback):
        return