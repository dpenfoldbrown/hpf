#!/usr/bin/env python

import os
import sys
import getopt
from hpf.runtime import runtime, Flag, IntegerOption, FileOption
from hpf.utilities import consume
from hpf.processing import processor, SYNCHRONOUS
from hpf.amnh.tree import TreeSplitter, UnrootedShortestPath
from Bio.Nexus.Nexus import Nexus
from Bio.Nexus.Trees import Tree
from Bio import AlignIO
from Bio.Align.Generic import Alignment
from StringIO import StringIO
from itertools import repeat
NEXUS = "nexus"
MAX_SIZE = "leaves"
DIRECTORY = "directory"

def _split(args):
    return split(*args)

def split(tree_file, size, nexus=False, dir=None):
    print file,size
    if nexus:
        tree = Nexus(tree_file).trees[0]
        tree2 = Nexus(tree_file).trees[0]
    else:
        with open(tree_file) as handle:
            tree_str = handle.read()
            tree = Tree(tree_str)
            tree2 = Tree(tree_str)
#    with open(align_file) as handle:
#        alignment = AlignIO.read(handle, "phylip")
    splitter = TreeSplitter(tree,max_size=size,annotater=UnrootedShortestPath)
    subs = list(splitter.subtrees())
    runtime().debug("Found",len(subs),subs)
    dir = dir if dir else os.path.dirname(tree_file)
    
    for i,tree in enumerate(subs):
        nodes = [tree.node(node) for node in tree.all_ids()]
        taxa = set([node.data.taxon for node in nodes if node.data.taxon != None])
        for terminal in tree2.get_terminals():
            node = tree2.node(terminal)
            if node.data.taxon in taxa:
                node.data.taxon = "%i-"%i + node.data.taxon
#        sub_taxa = tree.get_taxa()
#        sub_alignment = Alignment(alphabet=alignment._alphabet)
#        sub_alignment._records = [r for r in alignment._records if r.id in sub_taxa]
#        assert len(sub_taxa)==len(sub_alignment._records)
##        align_out = "%s.%i" % (os.path.join(dir,os.path.basename(align_file)),i)
#        with open(align_out,"w") as handle:
#            AlignIO.write([sub_alignment], handle, "phylip")
#        from hpf.phylip import interleave
#        interleave(align_out)
        with open("%s.%i" % (os.path.join(dir,os.path.basename(tree_file)),i),"w") as handle:
            print >>handle, tree.to_string(plain_newick=True,branchlengths_only=False)+";"
    with open("%s.annotated" % os.path.join(dir,os.path.basename(tree_file)),"w") as handle:
        print >>handle, tree2.to_string(plain_newick=True,branchlengths_only=False)+";"

def main(*args):
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS))
    runtime().debug("Using processor",pool)
    pool.make_tasks(lambda: zip(args,
#                                [args[i] for i in range(0,len(args),2)],
#                                [args[i] for i in range(1,len(args),2)],
                                repeat(runtime().opt(MAX_SIZE),len(args)),
                                repeat(runtime().opt(NEXUS),len(args)),
                                repeat(runtime().opt(DIRECTORY),len(args))
                                ))
    consume(pool.run(_split))
    
def _do(argv):
    r = runtime()
    r.description("""
    template.py [-options] args
    template script.
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing"))
    r.add_option(Flag(NEXUS, "n", description="The input files are in NEXUS format."))
    r.add_option(IntegerOption(MAX_SIZE, "l", description="Maximum number of leaves in the trees. Default:40", default=40))
    r.add_option(FileOption(DIRECTORY, "o", description="Directory to export to. Default to dirname(tree_file)", default=None))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
