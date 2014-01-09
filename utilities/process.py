'''
Created on Apr 28, 2010

@author: Patrick
'''
import subprocess
import time
class Process(object):
    """
    Context manager for handling subprocesses.
    """

    def __init__(self,**kwargs):
        self.kwargs = kwargs
    
    def _args(self):
        assert False, "Abstract method"
    
    def check_call(self):
        with self as popen:
            popen.communicate()
            assert popen.poll() == 0

    def call(self):
        with self as popen:
            popen.communicate()
            return popen.poll()
    
    def __enter__(self):
        self.popen = subprocess.Popen(self._args(),**self.kwargs)
        return self.popen
    
    def __exit__(self, type, value, traceback):
        """
        Wait until the process is finished.
        """
        if self.popen.poll() == None:
            self.popen.communicate()

import sys
import os
from hpf.runtime import runtime, INDENT, DEBUG

class Task(object):
    
    PRE=0
    RUN=1
    POST=2
    
    def __init__(self, force=False, runtime=None):
        """
        @param force: Force this task to be run, regardless of status. 
        @param run: A runtime object for overriding logs/debug/etc.
        """
        self.force = force
        self.result = None
        from copy import copy
        self._runtime = runtime
    
    def _classname(self):
        return self.__class__.__name__
    
    def run(self, force=None):
        """
        Run this task.  Will skip _pre,_do,_post processing if result exists
        and force is not True.
        You may override the object's force attribute directly within this method.
        """
        self.runtime().debug(self._classname(),"run()")
        force = force if force!=None else self.force
        if self.result and not force:
            self.runtime().debug("Already run.",INDENT=1)
            return self.result

        try:
            # Pre-Processing
            try:
                Task._pre(self)
                self._pre()
            except Exception as e:
                if self._except(Task.PRE, e):
                    raise
            # Standard processing
            try:
                self.runtime().indent()
                if self._skip() and not force:
                    Task._passed(self)
                    self.result = self._passed()
                else:
                    Task._do(self)
                    self.result = self._do()
            except Exception as e:
                if self._except(Task.RUN, e):
                    raise
            finally:
                self.runtime().unindent()
        finally:
            # Post processing
            try:
                Task._post(self)
                self._post()
            except Exception as e:
                if self._except(Task.POST, e):
                    raise
        
        return self
    
    def _do(self):
        """The driver method."""
        self.runtime().debug("Running",time.strftime("%H:%M:%S %a %m/%d/%y", time.localtime()),DEBUG=2)
        pass
    def _skip(self):
        """Whether this task is to be skipped."""
        pass
    def _passed(self):
        """Method to be called if the job is passed on and skipped."""
        self.runtime().debug("Skipping processing",time.strftime("%H:%M:%S %a %m/%d/%y", time.localtime()), DEBUG=2)
    def _post(self):
        """Post-processing method."""
        self.runtime().debug("Post-Process",time.strftime("%H:%M:%S %a %m/%d/%y", time.localtime()),DEBUG=2)
    def _pre(self):
        """Pre-processing method."""
        self.runtime().debug("Pre-Process",time.strftime("%H:%M:%S %a %m/%d/%y", time.localtime()),DEBUG=2)
    def _except(self,id,exception=None):
        """
        Method called if pre, post, or processing hits exception.
        @return: True/False continue to raise this exception.
        """
        self.runtime().debug(id,"Exception",exception,time.strftime("%H:%M:%S %a %m/%d/%y", time.localtime()))
        return True
    
    def runtime(self):
        """
        @return: Task specific runtime if available, otherwise global.
        """
        return self._runtime if self._runtime else runtime()


class TaskTree(Task):
    """
    A task that takes input arguments.  If the arguments in input
    are Tasks, they will be run by pre-processing.
    """
    
    def __init__(self,*input,**kwargs):
        super(TaskTree,self).__init__(**kwargs)
        self.input = []+list(input)
    
    def append_input(self, input):
        self.input.append(input)
        return input
    
    def _pre(self):
        """Process all input tasks."""
        self.input_results = []
        self.runtime().debug("Pre-Process required tasks")
        for task in self.input:
            if isinstance(task, Task):
                self.runtime().debug("Running required task",task._classname(),INDENT=1)
                result = task.run()
            else:
                result = task
            self.input_results.append(result)
        super(TaskTree,self)._pre()
        
        
class TaskHandler(TaskTree):
    """Handler class that can be used to run a chain or tree of tasks"""

    def __init__(self,dir=None,run=None,*args,**kwargs):
        """
        @param dir: Directory to run this task in. 
        """
        super(TaskHandler,self).__init__(*args,**kwargs)
        self.dir = dir
        
    def _pre(self):
        """Change directories before running anything."""
        if self.dir:
            self._cwd = os.getcwd()
            os.chdir(self.dir)
            self.runtime().debug("Changing dir:",self.dir)
        super(TaskHandler,self)._pre()
            
    def _post(self):
        """Run all post-processing and switch back to original directory."""
        super(TaskHandler,self)._post()
        if self.dir:
            os.chdir(self._cwd)
            self.runtime().debug("Reverting dir:",self.dir)

    def _do(self):
        self.runtime().debug("Finished")
        return self

    def add_task(self,*task):
        """
        @return: task, for chaining arguments
        """
        super(TaskHandler,self).append_input(*task)
        return task if len(task)!=1 else task[0]

class OutputTask(TaskTree):
    """Defines a simple task with input and output files that will skip execution if the file already exists."""

    def __init__(self,*output,**kwargs):
        super(OutputTask,self).__init__(**kwargs)
        self.output = list(output)
        if len(self.output) == 0:
            self.output = []
        elif len(self.output) == 1:
            self.output = self.output[0]
    
    def append_output(self, *output):
        if not hasattr(self.output,'__iter__') and len(output)>0:
            self.output = [self.output]
        self.output += output
        return self
            
    def _skip(self):
        output = self.output
        if not hasattr(output,'__iter__'):
            output = [output]
        self.runtime().debug("Required output files:")
        self.runtime().indent()
        for file in output:
            assert file!=None, "File: %s" % str(file)
            exists = os.path.exists(file)
            self.runtime().debug(file,"exists" if exists else "doesn't exist")
        self.runtime().unindent()
        return all([os.path.exists(f) for f in output])

    def _except(self, id, exception=None):
        if id==Task.RUN:
            self.runtime().debug("Processing failed, removing files:")
            self.runtime().indent()
            output = self.output
            if not hasattr(output,'__iter__'):
                output = [output]
            for file in [file for file in output if isinstance(file, basestring)]:
                try:
                    if os.path.exists(file):
                        os.remove(file)
                        self.runtime().debug(file)
                except:
                    self.runtime().printerr("Cannot remove file ",file)
                    raise
            self.runtime().unindent()
        return super(OutputTask,self)._except(id)

