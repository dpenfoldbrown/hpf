from optparse import OptionParser, OptionGroup
import sys
import os
from hpf.runtime import runtime

def consume(generator):
    """
    Consume a generator destroying items.
    @return: The number of items generated.
    """
    count = 0
    for r in generator:
        count +=1
    return count

class OutputTask(object):
    """Defines a simple task with input and output files that will skip execution if the file already exists."""

    def __init__(self,*output):
        self.output = list(output)
        if len(self.output) == 0:
            self.output = None
        elif len(self.output) == 1:
            self.output = self.output[0]
            
    def run(self,force=False):
        output = self.output
        if not hasattr(output,'__iter__'):
            output = [output]
        which = [os.path.exists(f) for f in output]
        if not all(which) or force==True:
            try:
                r = self._do()
                return r if r != None else self.output
            except:
                #print >>sys.stderr, "failure. deleting ",self.output
                delete = [self.output] if isinstance(self.output, str) else self.output 
                for file in delete:
                    try:
                        if os.path.exists(file):
                            os.remove(file)
                    except:
                        print >>sys.stderr, "cannot remove file ",file
                        pass
                raise
        else:
            runtime().debug("Output exists, skipping",[f for i,f in enumerate(output) if which[i]==True])
            return self._exists()
    
    def _do(self):
        raise NotImplementedError()
    def _exists(self):
        return self.output


class QueueParser(object):

    def __init__(self, file):
        if isinstance(file,str):
            self._handle = open(file)
        else:
            self._handle = file
        self._queue = []

    def _peek(self):
        if len(self._queue) < 1:
            line = self._handle.readline().strip()
            self._queue.append(line)
        return self._queue[0]

    def _get(self):
        self._peek()
        return self._queue.pop(0)
    
    def __iter__(self):
        return iter(self.next, None)

    def _finished(self):
        return self._peek() == None or self._peek() == ""

class defaultdict(dict):
        def __init__(self, default_factory=None, *a, **kw):
            if (default_factory is not None and
                not hasattr(default_factory, '__call__')):
                raise TypeError('first argument must be callable')
            dict.__init__(self, *a, **kw)
            self.default_factory = default_factory
        def __getitem__(self, key):
            try:
                return dict.__getitem__(self, key)
            except KeyError:
                return self.__missing__(key)
        def __missing__(self, key):
            if self.default_factory is None:
                raise KeyError(key)
            self[key] = value = self.default_factory()
            return value
        def __reduce__(self):
            if self.default_factory is None:
                args = tuple()
            else:
                args = self.default_factory,
            return type(self), args, None, None, self.items()
        def copy(self):
            return self.__copy__()
        def __copy__(self):
            return type(self)(self.default_factory, self)
        def __deepcopy__(self, memo):
            import copy
            return type(self)(self.default_factory,
                              copy.deepcopy(self.items()))
        def __repr__(self):
            return 'defaultdict(%s, %s)' % (self.default_factory,
                                            dict.__repr__(self))

def find(pattern, dir=None):
    """Searches a directory path and yields regular expression matches."""
    dir = dir if dir else os.getcwd()
    import re
    regex = re.compile(pattern)

    for (path, dames, fnames) in os.walk(dir) :
        for fn in fnames:
            abs = os.path.abspath(join(path, fn))
            match = regex.search(abs)
            if match:
                yield abs
    
if __name__ == "__main__":
    runtime = Runtime(description="Example runtime.")
    parser = runtime.parser()
    parser.add_option("-c","--complete",action="store_true",default=False,help="allow this script to complete")
    runtime.parse_args()
    assert(runtime.options.complete)
    runtime.debug("Good, you ran it")
