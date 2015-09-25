#!/usr/bin/env python

# A script to read in a line-delim'ed list of n GO terms and to return a "parent matrix"
# representing those terms, where M(i,j) = 1 if i is a parent of j (0 otherwise)
# Outputs matrix to tab-delim'ed file

import argparse

def build_common_parent_matrix(go_terms, database):
    """Build a matrix representing common parentage between two GO terms.
    C(i,j) = 1 if GO term i and GO term j share a direct parent (parent 1 step away).
    Otherwise, C(i,j) = 0
    """
    from hpf.superfunc.go import parent_finder
    
    # Build dictionary of parents (avoid redundant querying). Key: GO term
    parent_dict = dict()
    for term in go_terms:
        parent_dict[term] = parent_finder(term, database)
    
    # Build matrix
    matrix = list()
    for i in go_terms:
        row = list()
        for j in go_terms:
            # If parents[i] and parents[j] intersect, i & j have common parents (matrix -> 1)
            if len(set(parent_dict[i]) & set(parent_dict[j])) != 0: row.append(1)
            else: row.append(0)
        matrix.append(row)
    return matrix

def build_parent_matrix(go_terms, database):
    """Build a "parent matrix" for given GO terms. Par(i,j) = 1 if GO term i is a 
    parent of GO term j. Otherwise Par(i,j) = 0
    """
    from hpf.superfunc.go import local_leaf_children_finder
    matrix = list()
    for i in go_terms:
        row = list()
        children = local_leaf_children_finder([i, ], src_db=database)
        for j in go_terms:
            if j in children: row.append(1)
            else: row.append(0)
        matrix.append(row)
    return matrix

def main():
    from hpf.superfunc.go import read_file, output_matrix, output_sparse_matrix
    print "Reading GO terms from infile {0}".format(args.go_file)
    go_terms = read_file(args.go_file)
    
    if args.type == "parent":
        print "Building parent matrix against database {0}".format(args.database)
        result_matrix = build_parent_matrix(go_terms, args.database)
    elif args.type == "common_parent":
        print "Building common parents matrix against database {0}".format(args.database)
        result_matrix = build_common_parent_matrix(go_terms, args.database)

    print "Writing parent matrix to file {0}".format(args.outfile)
    if args.store_sparse:
        output_sparse_matrix(result_matrix, args.outfile, val=1)
    else:
        output_matrix(result_matrix, args.outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='GO parent matrix creator')
    parser.add_argument("-g", "--go_file", action="store", dest='go_file', required=True, 
            help='A line-delimited file containing GO terms to form the matrix')
    parser.add_argument("-o", "--outfile", action="store", dest='outfile', required=True,
            help='The output file to contain the GO parent matrix')
    parser.add_argument("-d", "--database", action="store", dest='database', default='mygolite_012012',
            help='The database against which to query GO ancestral relationships')
    parser.add_argument("--type", action="store", dest='type', default="parent",
            choices=["parent", "common_parent"],
            help='The type of matrix to produce. parent is whether term i is a parent of j. common_parent is whether term i and j have parent(s) in common')
    parser.add_argument("-s", "--sparse", action="store_true", dest="store_sparse", default=False,
            help='A flag for storing matrices in sparse format')
    
    args = parser.parse_args() 
    main()

