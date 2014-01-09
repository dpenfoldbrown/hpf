# Quick graph functions

import networkx

def max_coverage(graph, node):
    if graph.node[node]['cscore'] != -1:
        return graph.node[node]['cscore']
    preds = graph.predecessors(node)
    if not preds:
        graph.node[node]['cscore'] = graph.node[node]['score']
        return graph.node[node]['cscore']
    max_pred = None
    for pnode in preds:
        pnode_score = max_coverage(graph, pnode)
        if not max_pred or max_pred < pnode_score:
            max_pred = pnode_score
    graph.node[node]['cscore'] = graph.node[node]['score'] + max_pred
    return graph.node[node]['cscore']

#class Path

def path_score(graph, path):
# Takes a list of nodes, returns sum of their scores
    score = 0
    for node in path:
        score += graph.node[node]['score']
    return score

def max_score_path(graph, node, path=[]):
# Returns the path with the highest summed node score
# NOTE: Could be much more efficient if 'path' was an object that kept it's sum score
# as it was being built (no re-iteration to score paths in max calc.)
    #print "\t Path: {0}".format(path)
    path = path + [node]
    preds = graph.predecessors(node)
    if not preds:
        return path
    max_path = None
    for pred_node in preds:
        pred_path = max_score_path(graph, pred_node, path=path)
        if not max_path or path_score(graph, max_path) < path_score(graph, pred_path):
            max_path = pred_path
    return max_path
    

def main():
    d = networkx.DiGraph()
    d.add_node(1, score=250, cscore=-1)
    d.add_node(2, score=360, cscore=-1)
    d.add_node(3, score=320, cscore=-1)
    d.add_node(4, score=560, cscore=-1)
    d.add_node(5, score=360, cscore=-1)
    d.add_node(6, score=240, cscore=-1)

    d.add_edge(1,3)
    d.add_edge(1,4)
    d.add_edge(1,5)
    d.add_edge(1,6)
    d.add_edge(1,3)
    d.add_edge(2,3)
    d.add_edge(2,5)
    d.add_edge(2,6)
    d.add_edge(3,5)
    d.add_edge(3,6)
    d.add_edge(4,6)

    for node in d.nodes():
        path = max_score_path(d, node)
        print "Max path from node {0}: ".format(node), path, 
        print "\tScore: ", path_score(d, path)
    
    #print "Max of 1: ", max_coverage(d, 1)
    #print "Max of 2: ", max_coverage(d, 2)
    #print "Max of 3: ", max_coverage(d, 3)
    #print "Max of 4: ", max_coverage(d, 4)
    #print "Max of 5: ", max_coverage(d, 5)
    #print "Max of 6: ", max_coverage(d, 6)


if __name__ == "__main__":
    main()
