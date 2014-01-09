import os
import errno
from hpf.runtime import debug

TARRED = "tar"
BZIPPED = "bz"
GZIPPED = "gz"

def compare(file1, file2):
    "Compares two files modification times and returns true if file1 is newer than file2"
    return os.path.getmtime(file1) > os.path.getmtime(file2);

def find(pattern, dir=os.getcwd()):
    """Searches a directory path and yields regular expression matches."""
    import re
    regex = re.compile(pattern)

    for (path, dames, fnames) in os.walk(dir) :
        for fn in fnames:
            debug(fn)
            abs = os.path.abspath(join(path, fn))
            match = regex.search(abs)
            if match:
                yield abs

def lock(path):
    pass
#    debug('acquiring lock %s' % path)
#    while os.path.exists(path):
#        time.sleep(5)
#    try:
#        os.mkdir(path)
#    except:
#        lock(path)
#    debug('\tqcquired lock %s' % path)

def unlock(path):
    pass
#    debug('releasing lock %s' % path)
#    try:
#        os.rmdir(path)
#    except:
#        debug("no lock to unlock")

def getDirectorySize(directory):
    dir_size = 0
    for (path, dirs, files) in os.walk(directory):
        for file in files:
            filename = os.path.join(path, file)
            dir_size += os.path.getsize(filename)
    return dir_size

def ensure(path, tryagain=5):
    """Attempts to make a path, fails on error.  tryagain will sleep for x seconds and try again (default=5, 0 means don't)."""
    if not exists(path):
        try:
            makeDirs(path)
        except:
            if tryagain > 0:
                try:
                    import time
                    time.sleep(tryagain)
                    makeDirs(path)
                except:
                    pass
    return existsOrFail(path)

def makeDirs(path):
    "Raises a PathException if the path cannot be created w/ os.error number as value"
    from hpf.utilities import debug
    debug("Ensuring: %s" % path)
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise PathException("Cannot make dir", e.errno)

def existsOrFail(path):
    if not os.path.exists(path):
        from hpf.utilities import error
        error( "Path doesn't exist: %s" % path, path)
    return path

def exists(path):
    return os.path.exists(path)

def split(path):
    "Returns a touple (directory, file) of a path."
    return os.path.split(path);

def getFile(path):
    "Splits a path and returns the file"
    head, tail = os.path.split(path)
    return tail

def getDirectory(file):
    "Splits a path and returns the directory of the file"
    head, tail = os.path.split(file)
    return head

def getExtension(path):
    "Returns the file extension in a touple (ext, root)."
    root, ext = os.path.splitext(path)
    return ext, root

def removerf(path):
    from hpf.utilities import system
    system("rm -rf %s" % path)
#    "Dangerous! Recursively forces deletion of all files!"
#    # Delete everything reachable from the directory named in 'top',
#    # assuming there are no symbolic links.
#    # CAUTION:  This is dangerous!  For example, if top == '/', it
#    # could delete all your disk files.
#    for root, dirs, files in os.walk(top, topdown=False):
#        for name in files:
#            os.remove(os.path.join(root, name))
#        for name in dirs:
#            os.rmdir(os.path.join(root, name))


def rename(src, dest):
    os.rename(src, dest)
    return dest

def size(path):
    try:
        return os.path.getsize(path)
    except:
        return None

def join(*paths):
    return os.path.join(*paths)

class PathException(Exception):
    def __init__(self, message, value=None):
        self.message = message
        if value == None:
            value = message
        self.value = value
    def __str__(self):
        return repr(self.message)

def zipSuffix(path):
    "Return a touple, (zipType, suffix) or (False, None)"
    import re
    match = re.compile("\.(bz|gz)$", re.IGNORECASE).search(path)
    if match :
        suffix = match.string[match.start():match.end()]
        from hpf.utilities import conditional
        type = conditional(suffix.lower() == BZIPPED, BZIPPED, GZIPPED)
        return type, suffix
    else:
        return None, None

def isZipped(path):
    return zipSuffix(path)[0] != None

def isTarred(path):
    return tarSuffix(path)[0] != None

def tarSuffix(path):
    "Returns true if this path is tarred.  Returns a touple: (tarType, suffix)"
    import re
    match = re.compile("\.(tar\.|t)(bz|gz)$", re.IGNORECASE).search(path)
    if match :
        suffix = match.string[match.start():match.end()]
        z = re.compile("bz$|gz$", re.IGNORECASE).search(suffix)
        if z :
            from hpf.utilities import conditional
            type = conditional(z.string[z.start():z.end()].lower() == BZIPPED, BZIPPED, GZIPPED )
	else :
	    type = TARRED
        return type, suffix
    else:
        return None, None

