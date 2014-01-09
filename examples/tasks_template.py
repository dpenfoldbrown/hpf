#!/usr/bin/env python

# A test to work out the tasks/processing pool/make_tasks structure, and how it can be used (how
#   it works with different types of inputs, etc).
# dpb, 1/17/2011

from hpf.utilities import consume
from hpf.processing import processor

# Imports for simulation only
from glob import glob
from os.path import expanduser
from subprocess import call


def process(input):
    print "process:: Input: {0}".format(input)
    
    # Do some processing of input task.
    call(["head", "-n 1", input])
    

# Note that tasks function can take anything, as it is just called on the same args given to make_tasks.
def tasks(input_str):
    print "tasks:: Input: {0}".format(input_str)

    # Parse input string into list or iterable rep'ing individual elements.
    file_list = glob(expanduser(input_str))
    print "tasks:: Files: ", file_list
    
    # Return list/iterable/sequence/variable rep'ing "tasks" to give to process function.
    return file_list
    

def main():
    # Simulated input string (dir containing fasta files).
    input_str="/Users/dpb/bonneau-dev/sandbox/frag-input/oi_test/*.fasta"

    # Create a processing pool (via hpf.processing) and make tasks to serve to the process function.
    pool = processor()
    
    # Input to make_tasks can be any number of arguments, as they are passed directly to the given function.
    pool.make_tasks(tasks, input_str)
    consume(pool.run(process))


if __name__ == "__main__":
    main()    
    

# NOTE:
# The make_tasks function from hpf.processing takes a function and a variable number of inputs.
#   EG. make_tasks(func, *args)
# make_tasks will call the given function on the given input set, convert the returned value to a list,
# and store the list in self._tasks:
#   self._tasks = list(func(*args))
# So, the function given to make_tasks should return ANYTHING that can be converted into a list via
# python.list(), where the list elements represent "tasks" -> eg., input into the processor function
# that you pass to processor_pool.run(...).
#
# EG:
# pool.make_tasks(tasks, *files)
# pool.run(process_files)
# Where tasks returns an iterable over a number of files, or a list of files, or something, and
# process_files takes a filename as input.
    
# NOTE:
# Not using custom runtime functionality, because generally I think it is unnecessary. 
    
# NOTE:
# *args in a function definition (eg. def func(*args)) indicates func takes any number of non-keyworded
#   arguments.
# Calling a function with *list (eg. func2(*list)) calls the function with arguments being all items
#   in the list. EG:
#   list = [1,2,3]
#   func2(*list) is equivalent to func2(list[0], list[1], list[2])
