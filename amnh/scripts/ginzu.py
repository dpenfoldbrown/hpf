#!/usr/bin/env python

import os
import sys
import getopt
from numpy import mean
from collections import defaultdict
from hpf.runtime import runtime, Flag, debug
from hpf.utilities import consume
from hpf.processing import processor, SYNCHRONOUS

import pylab
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import matplotlib.font_manager as font
import networkx as nx
from itertools import izip
from numpy import arange

from families import annotate_records, domain_types, DOMAIN_TYPES
colors = defaultdict(lambda: '#AEBC21',
                    {'pfam':'#757116','heuristic':'#AEBC21','msa':'#D9DB56','fold_recognition':'#00477F','psiblast':'#4C88BE'})
        

def normalize(yield_dict):
    """Normalize the average domain coverage values."""
    y_dict = {}
    for key in yield_dict.keys():
        y_dict[key] = yield_dict[key]/sum(yield_dict.values())
    return y_dict

def family(fasta):
    """
    Annotate a fasta with family coverage
    @return: (family_name, coverage_dict)
    """
    base = os.path.basename(fasta)
    family_name = ".".join(base.split(".")[:-2])

    records = annotate_records(fasta)
    domains = domain_types(records)
    averages = dict()
    for type in DOMAIN_TYPES:
        averages[type] = mean(domains[type]) if len(domains[type])>0 else 0.0 
    return (family_name,normalize(averages))

def plot(families, filename="ginzu.svg"):
    #pylab.figure(figsize=(7,7))
    families.sort(cmp=lambda x,y: cmp(x[1]['psiblast'], y[1]['psiblast']))
    
    fp = font.FontProperties(size="x-small")
    ax = pylab.axes([0.3, 0.0, .6, .7])
    ind = arange(len(families))        
    width = .35
    plots = {}
    domain_dict = {}
    for type in DOMAIN_TYPES:
        domain_dict[type] = []
    
    names = []
    for name,domains in families:
        names.append(name)
        for type in DOMAIN_TYPES:
            domain_dict[type].append(domains[type])

    sum_bottom = [0 for i in arange(len(families))]
    debug("length of sum_bottom:", len(sum_bottom))
    for dkey in DOMAIN_TYPES:
        #plots[dkey] = pylab.bar(ind, domain_dict[dkey], color=self.color[dkey],bottom=sum_bottom)
        plots[dkey] = pylab.barh(ind, domain_dict[dkey], color=colors[dkey],left=sum_bottom)
        sum_bottom = map(sum, zip(sum_bottom,domain_dict[dkey]))
        
    #pylab.xticks(ind+width/2, organisms_dict.keys(), rotation=45)
    pylab.yticks(ind+width/2, names, size='xx-small')
    pylab.xticks()
    pylab.title("Ginzu domain frequencies")
    
    pylab.legend( [plots[key][0] for key in DOMAIN_TYPES], DOMAIN_TYPES, prop=fp, markerscale=.5, loc=(.85,.85))
    #kdrew: a little filename manipulation
    #f_parts = self.yield_plot_filename.rpartition('.')
    #filename = f_parts[0] + "_" + self.experiment_key + "." + f_parts[2]
    pylab.savefig(filename)

def main(*args):
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS))
    runtime().debug("Using processor",pool)
    pool.make_tasks(lambda: args)
    families = [f for f in pool.run(family) if not isinstance(f,Exception)]
    plot(families)

    
def _do(argv):
    r = runtime()
    r.description("""
    ginzu.py [-options] args
    Create a large graphic showing ginzu info for all protein families.
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing"))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
