'''
Created on Nov 23, 2009

@author: patrick
'''
import sys
import os
import getopt
from Queue import Queue

DEBUG = "debug"
SILENT = "silent"
HELP = "help"
INDENT = "indent"

def runtime():
    """
    @return: Runtime singleton
    """
    if not Runtime._instance:
        Runtime._instance = Runtime()
    return Runtime._instance

def debug(*args,**kwargs):
    runtime().debug(*args,**kwargs)

def println(*args):
    runtime().println(*args)


class Option(object):
    
    def __init__(self,name,short_option, expects=True, description=None, default=None):
        self.name = name
        self.short_option = short_option
        self.expects = expects
        self.description = description
        self.value = default
        assert len(self.short_option) == 1

    def set_value(self, value=None):
        self.value = value
        return self
        
    def validate(self,value):
        self.set_value(value)
        return self.value
    
    def short(self):
        return self.short_option+":" if self.expects else self.short_option

class Flag(Option):
    
    def __init__(self,name,short_option,description,**kwargs):
        Option.__init__(self,name,short_option,expects=False,**kwargs)
        
    def validate(self, value):
        assert value==None or value=="", "value is "+str(value)+"."
        return Option.validate(self, True)
        
class FileOption(Option):

    def validate(self, value):
        assert os.path.exists(value)
        return Option.validate(self, value)
    
class IntegerOption(Option):
    
    def validate(self,value):
        return Option.validate(self, int(value))

class Help(Option):
    
    def __init__(self):
        Option.__init__(self, HELP,"?", expects=False, description="Print out usage and help information.")
        
    
    def validate(self, value):
        r = runtime()
        print r.usage()
        sys.exit(1)

class Runtime(object):
    '''
    Class for managing debugging/printing, options, etc.
    '''
    
    _instance = None
    
    def __init__(self):
        self._short = {}
        self._long = {}
        self.add_option(IntegerOption(DEBUG, 'd', description="Debug level, requires int value.").set_value(0))
        self.add_option(Flag(SILENT, 'q', description="Silent mode.").set_value(False))
        self.add_option(Help())
        #self.add_option(option)
        self.dir_stack = list()
        self._log = None
        self._elog = None
        self._indent = 0
    
    def pushd(self, dir):
        cur_dir = os.getcwd()
        self.dir_stack.append(cur_dir)
        os.chdir(dir)
        return cur_dir
    
    def popd(self):
        dir = self.dir_stack.pop()
        os.chdir(dir)
        return dir
    
    def parse_options(self, argv):
        short_str = self._short_str()
        long_list = self._long_list()
        opts, args = getopt.getopt(argv, short_str, long_list)
        for o,a in opts:
            #print o,a
            option = self._short[o] if self._short.has_key(o) else self._long[o]
            option.validate(a)
        self._args = args
        return args
            
    def add_option(self, option):
        self._long['--'+option.name] = option
        self._short['-'+option.short_option] = option
        
    def description(self, value):
        self._description = value

    def usage(self):
        out = ""
        for line in self._description.split("\n"):
            if line=="":
                continue
            out += line.strip()+"\n"
        for key in sorted(self._short.keys()):
            option = self._short[key]
            out+="  -%s, --%s%s\t%s\n" % (option.short_option,option.name, "=VALUE" if option.expects else "", option.description if option.description != None else "")
        return out

    def _short_str(self):
        return "".join([self._short[key].short() for key in sorted(self._short.keys())])
    
    def _long_list(self):
        return [self._short[key].name for key in sorted(self._short.keys())]

    def set_log(self, stdout, stderr=None):
        self._log = stdout
        self._elog = stdout if stderr==None else stderr
    
    def log(self, *args,**kwargs):
        if self._log:
            print >>self._log, self._tabs(**kwargs)+" ".join([str(a) for a in args])
        
    def elog(self, *args,**kwargs):
        if self._elog:
            print >>self._elog, self._tabs(**kwargs)+" ".join([str(a) for a in args])
        

    def _opt_dump(self):
        return [(name, self.opt(name)) for name in self._long]
    
    def opt(self,name):
        """
        Retrieve an option value.
        """
        return self._long["--"+name].value

    def set_debug(self,level=1):
        """
        Set the global debug level.
        """
        option = self._long["--"+DEBUG]
        tmp = option.value
        option.set_value(level)
        return tmp

    def debug(self,*args,**kwargs):
        """
        Prints the arguments if the runtime debug level is >= 0|kwargs[DEBUG]
        @param kwargs: kwargs[DEBUG] can be set for this statement
        """
        if kwargs.has_key(DEBUG):
            if self.opt(DEBUG)>=kwargs[DEBUG]:
                self.println(*args,**kwargs)
        else:
            if self.opt(DEBUG)>0:
                self.println(*args,**kwargs)
    
    def _tabs(self,**kwargs):
        extra = kwargs[INDENT] if kwargs.has_key(INDENT) else 0
        return "".join(["\t" for i in xrange(self._indent+extra)])
    
    def indent(self, indent=None):
        """
        Indent printing by tabs.
        @type indent: integer number of tabs to indent.  
        """
        indent = indent if indent else 1
        self._indent += indent
        return self._indent - indent
    
    def unindent(self, indent=None):
        indent = indent if indent else 1
        self._indent -= indent
        return self._indent + indent
    
        
    def println(self,*args,**kwargs):
        """
        Prints the arguments if the runtime silent option is not set.
        """
        if not self.opt(SILENT):
            print self._tabs(**kwargs)+" ".join([str(a) for a in args])
        if self._log:
            self.log(*args,**kwargs)

    def printerr(self,*args,**kwargs):
        """
        Prints the arguments if the runtime silent option is not set.
        """
        if not self.opt(SILENT):
            print >>sys.stderr, self._tabs(**kwargs)+" ".join([str(a) for a in args])
        if self._elog:
            self.elog(*args,**kwargs)

