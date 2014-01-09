#!/usr/bin/env python

# A script to create annotation matrices for a set of MGI IDs and a list of
# GO terms.
# Matrix M(i, j) = 1 iff MGI ID i is annotated with GO term j in the given
# or set mygolite database. Can choose to include IEA annotations or not.

import argparse

def build_annotation_matrix(mgi_ids, go_terms, database):
    from hpf.superfunc.go import mgi_go_annotations
    matrix = list()
    for mgi in mgi_ids:
        row = list()
        mgi_terms = mgi_go_annotations(mgi, database=database, no_iea=True)
        for go_term in go_terms:
            if go_term in mgi_terms:
                row.append(1)
            else:
                row.append(0)
        matrix.append(row)
    return matrix


def main():
    from hpf.superfunc.go import read_file, output_matrix, output_sparse_matrix
    print "Reading MGI IDs ({0}) and GO terms ({1}) for matrix creation".format(args.mgi_file, args.go_file)
    mgi_ids  = read_file(args.mgi_file)
    go_terms = read_file(args.go_file)

    print "Creating annotation matrix against database {0}".format(args.database)
    annotation_matrix = build_annotation_matrix(mgi_ids, go_terms, args.database)

    print "Writing annotation matrix to file {0}".format(args.outfile)
    if args.store_sparse:
        output_sparse_matrix(annotation_matrix, args.outfile, val=1)
    else:
        output_matrix(annotation_matrix, args.outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MGI GO annotation matrix creator')
    parser.add_argument("-m", "--mgi_file", action="store", dest='mgi_file', required=True, 
            help='A line-delimited file containing MGI IDs to form the rows of the matrix')
    parser.add_argument("-g", "--go_file", action="store", dest='go_file', required=True, 
            help='A line-delimited file containing GO terms to form the columns of the matrix')
    parser.add_argument("-d", "--database", action="store", dest='database', default='mygolite_012012',
            help='The database against which to query MGI IDs to find GO annotations')
    parser.add_argument("-s", "--sparse", action="store_true", dest="store_sparse", default=False,
            help='A flag for storing matrices in sparse format')
    parser.add_argument("-o", "--outfile", action="store", dest='outfile', required=True,
            help='The output file containing MGI->GO term annotation matrix')
    args = parser.parse_args()
    main()
