import sys
import os
import cPickle
import datetime
import traceback
from hpf.runtime import runtime
try:
    import multiprocessing
except:
    pass

TYPE = 'type'
SYNCHRONOUS = 'synchronous'
PROCESSORS = 'processors'
MULTIPROCESSOR = 'multiprocessor'
    
def processor(*args,**kwargs):
    """Create a processor object.  Attempts to auto-discover the mode based on
    `which qsub`.
    dpb: This is a bad idea. Parsing output of qsub is ghetto and unappreciated.
    Will recode as follows. Function will take a kwarg type=<TypeString>.
    If type not given, will default to using 'which qsub' as was previously done.
    Otherwise, will create object based on type string.
    """
    # Supported types:
    supported_types = ["PBS", "SGE", "Synchronous", "Multiprocessor"]
    
    def has_and_true(key):
        """Utility kwargs function"""
        return kwargs.has_key(key) and kwargs[key]==True
    def find_num_processors():
        """Utility function to find number of processors for multiprocessing"""
        from numpy.distutils import cpuinfo
        if isinstance(cpuinfo.cpu.info,list):
            num_procs = len(cpuinfo.cpu.info)
        else:
            num_procs = int(cpuinfo.cpu.info['sysctl_hw']['hw.availcpu'])
        runtime().debug("Auto-discovered %i processors" % num_procs)
        return num_procs
    
    type = kwargs.pop(TYPE, False)
    # If type not given: use old 'which qsub' system. Otherwise use new type system
    if not type:
        # Check for synchronous
        synchronous = kwargs.pop(SYNCHRONOUS, False)
        if synchronous:
            return MapProcessor(*args,**kwargs)
        # If no synchronous, run which qsub
        import subprocess
        qsub = subprocess.Popen("which qsub", shell=True, stdout=subprocess.PIPE).communicate()[0].strip()
        if qsub == "" or has_and_true(MULTIPROCESSOR):
            proc = MultiProcessor
            if not kwargs.has_key(PROCESSORS) or (kwargs.has_key(PROCESSORS) and kwargs[PROCESSORS]==None):
                kwargs[PROCESSORS] = find_num_processors()
        elif "grid" in qsub:
            proc = SGEArrayProcessor
        elif "pbs" in qsub or "torque" in qsub:
            proc = PBSArrayProcessor
    else:
        if type == "PBS":
            proc = PBSArrayProcessor
        elif type == "SGE":
            proc = SGEArrayProcessor
        elif type == "Synchronous":
            proc = MapProcessor
        elif type == "Multiprocessor":
            proc = MultiProcessor
            if not kwargs.has_key(PROCESSORS) or (kwargs.has_key(PROCESSORS) and kwargs[PROCESSORS]==None):
                kwargs[PROCESSORS] = find_num_processors()
        else:
            raise TypeError("Given type '{0}' not supported.\nSupported types are: {1}".format(type, supported_types))
    return proc(*args,**kwargs)
    
def array_processor(task_pickle,pbs=False):
    if pbs:
        return PBSArrayProcessor(task_pickle=task_pickle)
    else:
        return SGEArrayProcessor(task_pickle=task_pickle)
    
class Processor:
    def __init__(self):
        pass
    
    def make_tasks(self,func,*args):
        """Call func to return a list of tasks"""
        self._tasks = list(func(*args))
        return self._tasks
    
    def run(self, do, tasks, result=lambda x: x): 
        """Run the processor on functions do() and result() over the variable number of tasks."""
        abstract

class SGEArrayProcessor(Processor):
    """A processor for tasks based on the SGE array jobs.  Processes the task indexed at SGE_TASK_ID.  Will load a pickled array of tasks."""
    def __init__(self, task_pickle="tasks.pickle",**kwargs):
        Processor.__init__(self)
        self._task_pickle = task_pickle

    def make_tasks(self,func,*args):
        """Only produce tasks if they haven't been pickled yet."""
        self._tasks = None
        if self._task_pickle!=None and not os.path.exists(self._task_pickle):
            Processor.make_tasks(self, func, *args)
    
    def run(self, do, result=lambda *x: x):
        """
        Run differs here because this might represent an array task.
        If $SGE_TASK_ID is defined it will only process the appropriate task.
        Otherwise it will
        """
        if self._task_id():
            yield self._task(do,result)
        else:
            yield self._qsub()
    
    def _qsub(self):
        """
        Dump the tasks to a pickle and qsub an array job.
        """
        # dpb: catch error here - if this is run on a headnode (an init. run) when an old
        # pickle file exists in the cwd, make_tasks will be skipped and this method called,
        # previously resulting in a failure (no self._tasks defined). Now handled. Poor.
        if self._tasks == None:
            raise Exception("Error: Attempting to create pickle without setting " \
                "tasks. Delete old pickle file: {0} and run again.".format(self._task_pickle))
        
        cPickle.dump(self._tasks,open(self._task_pickle,"w"))
        
        cmd = " ".join(sys.argv)
        name = sys.argv[0]
        if self._tasks==None:
            self._tasks=[None]

        cmd = self._cmd(len(self._tasks),name,cmd)        
        print cmd
        
    def _task(self,do,result):
        """
        Load the tasks pickle and run a single task.
        """
        tasks = cPickle.load(open(self._task_pickle))
        if tasks==None:
            arg = self._task_id()
        elif isinstance(tasks, list):
            arg = tasks[self._task_id()-1]
        else:
            arg = tasks
        return result(do(arg))

    def _cmd(self,length,name,cmd):
        return "qsub -t 1-%i -j y -N %s -V -S `which python` -o /dev/null %s" % (length,name,cmd)

    def _task_id(self):
        try:
            return int(os.getenv("SGE_TASK_ID",None))
        except:
            return None
        
class PBSArrayProcessor(SGEArrayProcessor):
    def _task_id(self):
        try:
            return int(os.getenv("PBS_ARRAYID",None))
        except:
            return None

#    def _cmd(self,length,name,cmd):
#        cmd = """
#        for ((i=1;i<=%i;i++)); do
#           qsub -vTASK=$i -o /dev/null -V -j oe -m oe -q ser2 -S `which python` -N %s %s
#        done
#        """ % (length,name,cmd)
#        return cmd

def _map_wrap(args):
    return MapProcessor._wrap(*args)

def _wrap_functions(args):
    #print "WRAP", args
    instance, do, result, task = args
    return MultiProcessor._wrap_functions(instance, do, task, result)

def _return(args):
    return args

class MapProcessor(Processor):
    """Simple processor that uses the map function.  FOR DEBUGGING."""
    
    def __init__(self, raise_errors=True, **kwargs):
        Processor.__init__(self)
        self.exceptions = 0
        self.raise_errors = raise_errors
    
    def _wrap(self, do, result, task):
        try:
            r = do(task)
            r = result(r)
        except Exception, e:
            if self.raise_errors:
                raise
            r = e
            self.exceptions += 1
        return r
    
    def run(self, do, tasks=None, result=_return):
        if tasks == None:
            tasks = self._tasks
        if result == None:
            result = _return
        r = map(_map_wrap, [(self, do, result, task) for task in tasks])
        return r

class MultiProcessor(Processor):
    """Runs an array of tasks on a single machine making used of multiple processes.  It will keep track of progress for you if you like."""
    
    def __init__(self, processors=8, completed_out=sys.stdout, result_out=None, modulus=100, raise_errors=False, *args):
        Processor.__init__(self, *args)
        MultiProcessor.completed = multiprocessing.Value('i', 0)
        MultiProcessor.results = multiprocessing.Value('i', 0)
        MultiProcessor.exceptions = multiprocessing.Value('i',0)
        
        self.processors = processors

        # For printing result/completion info
        self.completed_out = completed_out
        self.result_out = result_out
        self.modulus = modulus
        self.raise_errors = raise_errors
        self.reset()
        self.pickle_file=None

    def reset(self):
        MultiProcessor.completed.value=0
        MultiProcessor.results.value=0
        MultiProcessor.exceptions.value=0

    def run(self, do, tasks=None, result=_return, batches=None):
        """Returns an unordered generator of results.  You can have the pools run in sets the size of 'batches' to avoid memory issues."""
        if tasks==None:
            tasks = self._tasks
        if len(tasks) < 1:
            raise Exception("No Tasks to Run", tasks, self.pickle_file)
        if result == None:
            result = _return
        self.total = len(tasks)
        self.start = datetime.datetime.now()
        if not batches:
            batches = len(tasks)
        print "Processor started with",len(tasks),"tasks at",self.start
        while(len(tasks)>0):
            batch = tasks[:batches]
            tasks = tasks[batches:]
            pool = multiprocessing.Pool(self.processors)
            for r in pool.map(_wrap_functions, [(self, do, result, task) for task in batch]):
                yield r
        #pool.join()
        pool.close()
        pool = None
        self.exceptions = MultiProcessor.exceptions.value

    def _wrap_functions(self, do, task, result):
        """Wrap the do/result functions so we can process and handle results in one pass.  By default, if result() is null will just return from do()"""
        try:
            try:
                r = do(task)
            finally:
                self._completed(task, None)
            r = result(r)
            self._result(task, r)
        except Exception, e:
            if self.raise_errors:
                traceback.print_exc(file=sys.stdout)
                raise
            else:
                traceback.print_exc(file=sys.stdout)
            r = e
            self.exceptions.value += 1
        return r

    def _completed(self, task, r):
        if not self.completed_out:
            return
        completed = MultiProcessor.completed
        completed.value +=1
        if completed.value % self.modulus == 0:
            now = datetime.datetime.now()
            delta = now-self.start
            s_per_task = float(delta.seconds)/float(completed.value)
            s_remaining = s_per_task*(float(self.total)-float(completed.value))
            print "tasks remaining:",(self.total-completed.value),int(s_remaining/60),"mins",r

    def _result(self, task, r):
        if not self.result_out:
            return
        results = MultiProcessor.results
        results.value +=1
        if results.value % self.modulus == 0:
            now = datetime.datetime.now()
            delta = now-self.start
            s_per_task = float(delta.seconds)/float(results.value)
            s_remaining = s_per_task*(float(self.total)-float(results.value))
            print "post remaining:",(self.total-results.value),int(s_remaining/60),"mins",r

# Static set up for typing
TYPES = {}
for i,type in enumerate([SGEArrayProcessor, PBSArrayProcessor, MultiProcessor]):
    classname = str().__class__.__name__
    TYPES[i] = type
    TYPES[classname] = type
    TYPES[type] = type
