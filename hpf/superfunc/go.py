#!/usr/bin/env python

# A set of GO related functions for superfuncck

def read_file(filename):
    """Reads and returns a list of values from line-delimited files"""
    values = list()
    with open(filename) as handle:
        for line in handle:
            values.append(line.rstrip())
    return values

def output_matrix(matrix, outfile):
    """Writes a matrix (list of lists) to given file. Tab-delimited"""
    handle = open(outfile, 'w')
    for row in matrix:
        for j in row:
            handle.write("{0}\t".format(j))
        handle.write("\n")
    handle.close()

def output_sparse_matrix(matrix, outfile, val=1):
    """Writes a matrix in sparse format. EG: <row>\t<col>\t<val>\n, where only entries with value='val'
    are written.
    """
    handle = open(outfile, 'w')
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] == val:
                handle.write("{0}\t{1}\t{2}\n".format(i, j, matrix[i][j]))
    handle.close()
        
def mgi_go_annotations(mgi_id, database='mygolite_012012', no_iea=True, host='handbanana.bio.nyu.edu', user='dpb', passwd='dpb_nyu'):
    """Takes an MGI ID (form MGI:1234567) and returns all GO annotations for that ID
    found in the given or default database. no_iea specifies whether IEA annotations
    are returned (by default, they are not returned).
    Other inputs are optional DB specifiers. Note that (as it is MGI), only mouse gene products will be queried.
    NOTE: Does not include negative associations (with association.is_not = 1)
    Returns empty list if no annotations are found
    """
    import MySQLdb
    db = MySQLdb.connect(host=host, user=user, passwd=passwd, db=database)
    cursor = db.cursor()

    query = "select t.acc " +\
            "from dbxref as d join gene_product as g on d.id=g.dbxref_id " +\
                "join association as a on g.id=a.gene_product_id " +\
                "join evidence as e on a.id=e.association_id " +\
                "join term as t on a.term_id=t.id " +\
            "where a.is_not=0 and d.xref_key='" + mgi_id + "'"
    if no_iea:
        query += " and e.code!='IEA'"

    count = cursor.execute(query)
    if count < 1: return list()
    results = cursor.fetchall()
    
    go_terms = set()
    for (go_term, ) in results:
        go_terms.add(go_term)
    return list(go_terms)

def parent_finder(go_term, database='mygolite_012012', host='handbanana.bio.nyu.edu', user='dpb', passwd='dpb_nyu'):
    """Takes a single go term string (form: GO:0001234). Returns the immediate parents
    (distance 1 UP the GO tree) of the term as found in the DB 'database'. Return is a 
    list of strings.
    """
    import MySQLdb
    db = MySQLdb.connect(host=host, user=user, passwd=passwd, db=database)
    cursor = db.cursor()

    parent_query = "select t.acc " +\
            "from graph_path as gp, term as t, term as t2 " +\
            "where gp.term1_id=t.id and gp.term2_id=t2.id and distance = 1 and t.acc != 'all' and t2.acc="
    
    parents = set()
    count = cursor.execute(parent_query + "'" + go_term + "'")
    if count < 1: return list()
    results = cursor.fetchall()
    for (parent, ) in results:
        parents.add(parent)
    return list(parents)


def local_leaf_finder(go_terms, src_db='mygolite_012012'):
# Takes a list of GO IDs as string (form: GO:00011222), returning the subset of given list containing
# the lowest GO terms on each branch represented by the given set (aka the local leaves)
    if not go_terms:
        print "No GO terms given"
        return None

    import MySQLdb
    db = MySQLdb.connect(host='handbanana.bio.nyu.edu', user='dpb', passwd='dpb_nyu', db=src_db)
    cursor = db.cursor()
    ancestors_query = "select t.acc " +\
            "from graph_path as gp, term as t, term as t2 " +\
            "where gp.term1_id=t.id and gp.term2_id=t2.id and distance != 0 and t.acc != 'all' and t2.acc="
    
    go_terms = set(go_terms)
    all_ancestors = set()
    for go in go_terms:
        count = cursor.execute(ancestors_query + "'" + go + "'")
        if count < 1: continue
        ancestors = cursor.fetchall()
        for (ancestor,) in ancestors:
            all_ancestors.add(ancestor)
    return list( go_terms.difference(all_ancestors) )

def local_leaf_children_finder(go_terms, src_db='mygolite_012012'):
# Takes a list of GO IDs as string (form: GO:00011222), returning all children of the list
    if not go_terms:
        print "No GO terms given"
        return None

    import MySQLdb
    db = MySQLdb.connect(host='handbanana.bio.nyu.edu', user='dpb', passwd='dpb_nyu', db=src_db)
    cursor = db.cursor()
    ancestors_query = "select t2.acc " +\
            "from graph_path as gp, term as t, term as t2 " +\
            "where gp.term1_id=t.id and gp.term2_id=t2.id and distance = 1 and t.acc != 'all' and t.acc="

    go_terms = set(go_terms)
    all_childz = set()
    for go in go_terms:
        count = cursor.execute(ancestors_query + "'" + go + "'")
        if count < 1: continue
        childz = cursor.fetchall()
        for (child,) in childz:
            all_childz.add(child)
    return list(all_childz)

