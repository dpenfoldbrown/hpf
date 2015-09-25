#!/usr/bin/env python

import os
import sys
import getopt
import MySQLdb
import networkx as nx
import subprocess
from hpf.bionet.export import XGMMLExporter
from hpf.bionet import Node,Edge,JNLP as DEFAULT_JNLP


DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
EXPERIMENT = "experiment"
NEIGHBORS = "neighbors"
SEQUENCES = "sequences"
FILENAME = "filename"
ORIGINAL = "original"
JNLP = "jnlp"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, EXPERIMENT:None, NEIGHBORS:False, SEQUENCES:None, FILENAME:None, ORIGINAL: False, JNLP: DEFAULT_JNLP}

__hpf = None
__err = None

HDDB_URL = "http://carl.bio.nyu.edu:5000/protein/%i"

def _graph(name,nodes,edges):
    graph = nx.MultiDiGraph(name=name)
    graph.add_nodes_from(nodes)
    for edge in edges:
        graph.add_edge(edge.get_source(),edge.get_target(),attr_dict=edge.attributes())
    return graph

def __key_gi(nodes):
    gis={}
    for n in nodes :
        if n.get("gis") != None and n.get("gis").strip() != "":
            #_debug("Name: ",n.get("name"))
            #_debug("PROTEIN: ",n.get("hddb"))
            for gi in n.get("gis").split(","):
                #_debug("GI: ",gi)
                try:
                    gis[int(gi)]=n
                except ValueError:
                    continue
    return gis

def _edges(nodes, sources, neighbors=False):
    """Grab interaction sources for a set of nodes.  If a node doesn't have
    a GI it is ignored"""
    _debug("Sources",sources)
    # Key gis for edge creation
    gis = __key_gi(nodes)
    assert(len(gis.keys())>0)
    # Get all the taxonomy_ids
    tax = set([n.get("taxonomy_id") for n in gis.values() if n.get("taxonomy_id") != None])
    assert(len(tax)==1)
    tax = list(tax)[0]
    
    # Query for interactions
    db = _err()
    cursor = db.cursor(MySQLdb.cursors.DictCursor)
    gi_list = ",".join([str(gi) for gi in gis.keys()])
    _debug(len(gis.keys())," GI's for edge searching")
    condition = "or" if neighbors else "and"
    for source in sources:
        _debug("Loading from ",source)
        if source=='kegg2':
            query = "select tablename from kegg2.species where taxid=%i" % tax
            cursor.execute(query)
            _debug(query)
            try:
                source = "%s.%s" % (source,cursor.fetchone()['tablename'])
                _debug("Kegg table",source)
            except:
                continue
        else:
            source = "%s.interactions" % source
        query = """
        select i1,i2, interactionType, direction from %s
        where i1 in (%s) %s i2 in (%s) 
        """ % (source,gi_list,condition,gi_list)#,tax,tax)
        #and taxid1=%i and taxid2=%i
        #_debug(query)
        cursor.execute(query)
        for row in cursor.fetchall():
            # Get Node by gi or mark it as neighbor with GI as name
            source = gis[row["i1"]] if gis.has_key(row["i1"]) else Node(row["i1"], {"original":0})
            target = gis[row["i2"]] if gis.has_key(row["i2"]) else Node(row["i2"], {"original":0})
            yield Edge((source, target), attr=row)
    cursor.close()
    
def query(experiment_key, table=None, sequences=None):
    """Query for node data from an individual experiment."""
    
    table_join = " join %s t on p.sequence_key=t.sequence_key " % table if table !=None else ""
    sequence_where = " and p.sequence_key in (%s) " % ",".join(sequences) if sequences != None else ""
    
    queries = []
    # Those with GI's
    # Must be done to get the specific GI
    query = """ 
    select 
    group_concat(DISTINCT a.gi SEPARATOR ',') as gis,
    b.size,b.pdb,b.rosetta,b.quality,b.sccs,b.ginzu,b.mf_value,b.mf_name,b.mf_acc,b.mf_pls_llr,b.mf_base_llr,b.coverage,b.id,
    p.experiment_key,p.sequence_key as parent_sequence_key,p.id as protein_key,
    e.bionet_taxonomy_id as taxonomy_id,
    m.acc as iea_acc, m.name as iea_name, ifnull(m.metric,0.0) as iea_llr,
    y.yrc_protein_key,
    ifnull(c.pmids,0) as pmid_count
    from experiment e 
    join protein p on e.id=p.experiment_key
    %s
    join sequenceAc a on p.sequence_key=a.sequence_key and a.taxonomy_id=e.bionet_taxonomy_id
    left outer join bionet_builder b on a.gi=b.gi
    left outer join protein_mf_min m on p.sequence_key=m.sequence_key
    left outer join yeastrcu y on p.sequence_key=y.sequence_key
    left outer join publications_count c on c.gi=a.gi
    where e.id=%i
    %s
    group by p.sequence_key;
    """  % (table_join, experiment_key, sequence_where)
    yield query
    # Those without GI's
    query = """ 
    select
    NULL as gis,
    b.size,b.pdb,b.rosetta,b.quality,b.sccs,b.ginzu,b.mf_value,b.mf_name,b.mf_acc,b.mf_pls_llr,b.mf_base_llr,b.coverage,b.id,
    p.experiment_key, p.sequence_key as parent_sequence_key,p.id as protein_key,
    e.bionet_taxonomy_id as taxonomy_id,
    m.acc as iea_acc, m.name as iea_name, ifnull(m.metric,0.0) as iea_llr,
    y.yrc_protein_key,
    ifnull(c.pmids,0) as pmid_count
    from experiment e 
    join protein p on e.id=p.experiment_key
    %s
    left outer join sequenceAc a on p.sequence_key=a.sequence_key and a.taxonomy_id=e.bionet_taxonomy_id
    left outer join bionet_builder b on p.sequence_key=b.parent_sequence_key
    left outer join protein_mf_min m on p.sequence_key=m.sequence_key
    left outer join yeastrcu y on p.sequence_key=y.sequence_key
    left outer join publications_count c on c.gi=a.gi
    where e.id=%i and a.gi is NULL
    %s
    group by p.sequence_key;
    """ % (table_join, experiment_key, sequence_where)
    yield query 
    

def _populate(edges,experiment_key):
    """Populate nodes for edges that were first neighbors but are not part of the initial node query.""" 
    db = _hpf()
    cursor = db.cursor(MySQLdb.cursors.DictCursor)
    cursor.execute("SET SESSION group_concat_max_len = 1000000;")
    
    gis = []
    for e in edges:
        for n in e:
            if n.get("original")==0:
                gis.append(n.get("name"))
    if len(gis) ==0:
        return
    gi_list = ",".join(set(gis))
    query = """
    select 
    group_concat(DISTINCT a.gi SEPARATOR ',') as gis,
    b.size,b.pdb,b.rosetta,b.quality,b.sccs,b.ginzu,b.mf_value,b.mf_name,b.mf_acc,b.mf_pls_llr,b.mf_base_llr,b.coverage,b.id,
    p.experiment_key,p.sequence_key as parent_sequence_key,p.id as protein_key,
    e.bionet_taxonomy_id as taxonomy_id,
    m.acc as iea_acc, m.name as iea_name, ifnull(m.metric,0.0) as iea_llr,
    y.yrc_protein_key,
    ifnull(c.pmids,0) as pmid_count
    from experiment e 
    join protein p on e.id=p.experiment_key
    join sequenceAc a on p.sequence_key=a.sequence_key and a.taxonomy_id=e.bionet_taxonomy_id
    left outer join bionet_builder b on p.sequence_key=b.parent_sequence_key
    left outer join protein_mf_min m on p.sequence_key=m.sequence_key
    left outer join yeastrcu y on p.sequence_key=y.sequence_key
    left outer join publications_count c on c.gi=a.gi
    where e.id=%i and a.gi in (%s)
    group by p.sequence_key;
    """ % (experiment_key,gi_list)
    nodes = list(_nodes([query]))
    _debug("Neighbor nodes loaded ",len(nodes))
    gis = __key_gi(nodes)
    # Make new edges with fully loaded nodes
    for edge in edges:
        source, target = edge
        try:
            if source.get("original") == 0:
                source = gis[int(source.get("name"))]
                source.set("original",0)
            if target.get("original") == 0:
                target = gis[int(target.get("name"))]
                target.set("original",0)
        except KeyError:
            continue
        yield Edge((source, target), attr=edge.attributes())

def _nodes(queries):
    """Perform queries to return nodes
    @return: Generator of Node objects built from attribute rows"""
    db = _hpf()
    cursor = db.cursor(MySQLdb.cursors.DictCursor)
    cursor.execute("SET SESSION group_concat_max_len = 1000000;")
    for query in queries:
        _debug(query)
        cursor.execute(query)
        for row in cursor.fetchall():
            row['psiblast'] = int(row['ginzu'].find('psiblast') != -1 if row['ginzu'] != None else False)
            row['fold_recognition'] = int(row['ginzu'].find('fold_recognition') != -1 if row['ginzu'] != None else False)
            row['hddb_url'] = HDDB_URL % row['protein_key']
            gis = row["gis"].split(",") if row["gis"] != None else None
            name = gis[0] if gis != None else "hddb:%i"%row["parent_sequence_key"]
            #if row['yrc_protein_key'] != None and is_int(row['yrc_protein_key']):
            #    row['url'] = "http://www.yeastrc.org/pdr/viewProtein.do?id=%i" % int(row['yrc_protein_key'])
            if gis!=None:
                row["url"] = "http://www.yeastrc.org/pdr/quickSearch.do?type=description&query=%s" % ",".join([gi for gi in gis if is_int(gi)])
            node = Node(name,attr=row)
            assert node.get("hddb") != None
            yield node
    cursor.close()

def is_int(s):
    try:
        int(s)
        return True
    except:
        return False

def _ename(experiment_key):
    with _hpf() as cursor:
        query = "select short_name from experiment where id=%i"%experiment_key
        _debug(query)
        cursor.execute(query)
        name = cursor.fetchone()[0]
        return name

def _hpf(db="hpf"):
    global __hpf
    if __hpf == None:
        __hpf = MySQLdb.connect(db=db,read_default_file="~/.my.cnf",host="127.0.0.1")
    return __hpf

def _err(db="bionetbuilder_info"):
    global __err
    if __err == None:
        __err = MySQLdb.connect(host="err.bio.nyu.edu",db=db,read_default_file="~/.my.cnf")
    return __err

def main():
    table=None
    seqs=None
    if opts[SEQUENCES]:
        if os.path.exists(opts[SEQUENCES]):
            with open(opts[SEQUENCES]) as f:
                seqs = [l.strip() for l in f]
        else:
            table = opts[SEQUENCES]
    
    experiments=[opts[EXPERIMENT]] if opts[EXPERIMENT] else _all_exp()
    for exp in experiments:
        _print("Experiment",exp)
        if (table or seqs) and opts[ORIGINAL]:
            try:
                each(exp,None,None,False)
            except Exception, e:
                _debug(e)
                pass
        try:
            each(exp,table,seqs,opts[NEIGHBORS])
        except Exception, e:
            _debug(e)
            pass

def _all_exp():
    with _hpf() as cursor:
        query = "select id from experiment where short_name is not NULL"
        _debug(query)
        cursor.execute(query)
        return [result[0] for result in cursor.fetchall()]

def each(experiment_key, table, seqs, neighbors=False):        
    nodes = list(_nodes(query(experiment_key,table,seqs)))
    for node in nodes:
        node.set("original",1)
    _print("Nodes",len(nodes))
    edges = list(_edges(nodes, ["bind1","biogrid3","dip3","intact3","kegg2"],neighbors=neighbors))
    _print("Edges",len(edges))
    # Populate nodes for neighbor edges
    if neighbors:
        neighbors = [edge for edge in edges if edge.get_source().get("original")==0 or edge.get_target().get("original")==0]
        original = [edge for edge in edges if edge.get_source().get("original")==1 and edge.get_target().get("original")==1]
        neighbors = list(_populate(neighbors,experiment_key))
        edges = neighbors+original
        nodes = set(nodes)
        for source,target in edges:
            nodes.add(source)
            nodes.add(target)
    
    _print("FINAL: ",len(nodes)," nodes ",len(edges)," edges")
    exp_name=_ename(experiment_key).replace(" ","_").replace(".","")
    graph = _graph(exp_name,nodes,edges)
    sub = "_sub" if table or seqs else ""
    file = opts[FILENAME] if opts[FILENAME] else "%s%s.xgmml" % (exp_name,sub)
    with XGMMLExporter(file) as exporter:
        exporter.write(graph,pretty_print=False)
    
    with open(os.devnull) as handle:
        jar_file = file+".jar"
        cmd = "jar cvf %s %s" % (jar_file, file)
        subprocess.check_call(cmd,shell=True,stdout=handle,stderr=handle)
    if opts[JNLP] != None:
        jnlp_file = file.replace("xgmml","jnlp")
        _print("JNLP",jnlp_file)
        with open(jnlp_file,"w") as w:
            with open(opts[JNLP]) as r:
                w.write(r.read()%(file,file))
    
    
def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "?dxqe:ns:f:oj:", [DEBUG,SILENT,UPLOAD,EXPERIMENT,NEIGHBORS,SEQUENCES,FILENAME,ORIGINAL,JNLP])
    global opts
    for o,a in _opts:
        if o in ('-?',HELP):
            _usage()
            sys.exit()
        if o in ('-d',DEBUG):
            opts[DEBUG] = True
        if o in ('-x', UPLOAD):
            opts[UPLOAD] = False
        if o in ('-q',SILENT):
            opts[SILENT] = True
        if o in ('-e',EXPERIMENT):
            opts[EXPERIMENT] = int(a)
        if o in ('-n',NEIGHBORS):
            opts[NEIGHBORS] = True
        if o in ('-s',SEQUENCES):
            opts[SEQUENCES] = a
        if o in ('-f',FILENAME):
            opts[FILENAME] = a
        if o in ('-o',ORIGINAL):
            opts[ORIGINAL] = True
        if o in ('-j',JNLP):
            assert os.path.exists(a)
            opts[JNLP] = a
            #_print(opts[JNLP])
    
    if opts[NEIGHBORS]:
        assert opts[SEQUENCES], "First neighbors don't make sense if we're doing the whole organism"
    return opts, args

def _debug(*args):
    if opts[DEBUG]:
        _print(*args)

def _print(*args):
    if not opts[SILENT]:
        print " ".join([str(a) for a in args])

def _usage():
    print "$ export_networks.py [-options] args"
    print "Export some bionetbuilder networks."
    print ""
    print "Options"
    print "-s sequence file or db table name. Usually some sub-selection of"
    print "    proteins like those without molecular function. 'sequence_key'"
    print "-o export the entire experiment as well as the sequence subset"
    print "-f filename to use"
    print "-j the jnlp file to use"
    print "-e build a specific experiment"
    print "-n attach first neighbors"


if __name__=="__main__":
    _do(sys.argv[1:])
