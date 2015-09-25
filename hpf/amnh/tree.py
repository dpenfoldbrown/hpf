'''
Created on Apr 19, 2010

@author: patrick
'''
from Bio.Nexus.Trees import Tree
from Bio.Nexus.Nodes import Node

class TreeAnnotater(object):
    
    def __init__(self,tree,max_size):
        self.tree = tree
        self.max_size = max_size
    
    def annotate(self):
        return self._annotate()
    
    def sorted(self):
        return self._sorted()
    
    def cost(self, root):
        return self._cost(root)
    
    def subtree(self, root):
        return self._subtree(root)

class RootedTerminals(TreeAnnotater):
    """
    Annotates nodes with distance to terminal information.
    """

    def _annotate(self):
        """
        Each node.data.terminals = # of terminals downstream from this node.
            leaf.data.terminals=1
        @return: self
        """
        tree = self.tree
        RootedTerminals._walk(tree)
        tree.count_terminals = lambda n: tree.node(n).data.terminals
        return self

    def _sorted(self):
        """
        @param kwargs: arguments passed to sorted() function. 
        @return: list nodes sorted by # of terminals
        """
        tree = self.tree
        terminals = lambda x: tree.node(x).data.terminals
        #ids = [id for id in tree.all_ids() if terminals(id)<=max_size]
        ids = tree.all_ids()
        compare = lambda x,y: cmp(terminals(x),terminals(y))
        return sorted(ids, cmp=compare)
        #return [id for id in sorted(tree.all_ids(), cmp=compare, **kwargs) if tree.node(id).data.terminals<=max_size]
    
    @staticmethod
    def _walk(tree,node=None):
        """
        Recursive depth-first walk of the tree annotating # of terminals.
        """
        if node is None:
            node=tree.root
        if tree.is_terminal(node):
            n = 1
        else:
            n = sum([RootedTerminals._walk(tree,succ) for succ in tree.node(node).succ])
        tree.node(node).data.terminals=n
        return n
    
    def _subtree(self,root):
        """
        Split a tree on a given node, pruning from the original tree.
        @param root: the node to use as the new root 
        @return: subtree rooted on this node.
        """
        sub = Tree(weight=self.tree.weight, 
                        rooted=self.tree.rooted, 
                        name=self.tree.name, 
                        data=self.tree.dataclass, 
                        max_support=self.tree.max_support)

        sub.node(sub.root).data = self.tree.node(root).data
        def _add(old_id,new_id):
            """
            Walk from this node, using the id from the old tree, and the id
            from the new tree to both load the data from the old tree and link 
            to the correct node in the new tree.
            """
            for old_succ in self.tree.node(old_id).succ:
                #print old_id,new_id
                to_add = Node(data = self.tree.node(old_succ).data)
                new_succ = sub.add(to_add,new_id)
                #print "\t",old_succ,new_succ
                _add(old_succ,new_succ)
        _add(root,sub.root)
        self.annotater.annotate(sub)
        unlink(self.tree,root)
        return sub
    
    def _cost(self, root):
        """
        The cost of a tree is defined as the #of terminals - max_size.
        """
        c = self.tree.node(root).data.terminals-self.max_size
        assert c<=0
        return abs(c)

class UnrootedShortestPath(TreeAnnotater):
    """
    Annotates nodes with shortest path to all other nodes
    """

    def _annotate(self):
        tree = self.tree
        self._g = self._graph()
        import networkx as nx
        length = nx.all_pairs_shortest_path_length(self._g)
        #print sorted(list(self._g.nodes()))
        #print sorted(list(self._g.edges()))
        for id in self._g:
            node = tree.node(id)
            node.data.length = length[id]

    def _shortest_path(self, x):
        """Shortest length to a terminal"""
        node = self.tree.node(x)
        terminals = self.tree.get_terminals()
        lengths = [self.tree.node(x).data.length[t] for t in terminals]
        return min(lengths)
    
    def _sorted(self):
        tree = self.tree
        terminals = tree.get_terminals()
        #internal = set(tree.all_ids())-set(terminals)
        ids = [id for id in self._g]
        return ids
        sorter = lambda x,y: cmp(self._shortest_path(x), self._shortest_path(y)) 
        return sorted(ids,cmp=sorter)
    
    def _graph(self):
        #print tree.unrooted
        import networkx as nx
        G=nx.Graph()
        tree = self.tree
        if len(tree.all_ids())==0:
            return G
        tree.unroot()
        def massage(branch):
            n1,n2,weight = branch
            return (n1,n2,weight)
        G.add_weighted_edges_from([massage(branch[0:3]) for branch in tree.unrooted])
        return G
        #>>> length=nx.all_pairs_shortest_path_length(G)
        #>>> print length[1][4]
        #>>> length[1]
        #{0: 1, 1: 0, 2: 1, 3: 2, 4: 3}

    def _cost(self, root):
        """
        Define the cost of a tree with the given root as the 
        distance to the Nth closest species.
        """
        g = self._subgraph(root)
        weight_sum = sum([edata['weight'] for u,v,edata in g.edges(data=True)])
        return weight_sum
        tree = self.tree
        closest_species = self._closest(root)
        return tree.node(root).data.length[closest_species[-1]]
    
    def _closest(self, id):
        """
        Return a list of the species ordered by distance from this node id.
        """
        node = self.tree.node(id)
        terminals = self.tree.get_terminals()
        sorter = lambda x,y: cmp(node.data.length[x], node.data.length[y]) 
        closest_species = sorted(terminals, cmp=sorter)
        return closest_species
    
    def _subgraph(self, root):
        import networkx as nx
        targets = self._closest(root)[0:self.max_size]
        # Get a directed tree
        g = nx.DiGraph()
        #t = nx.dfs_tree(self._g)
        # Add each path to target
        #print type(g)
        #print type(root)
        paths = nx.single_source_dijkstra_path(self._g,root)
        for target in targets:
            path = paths[target]
            edges_in_path = set()
            for i,n1 in enumerate(path[:-1]):
                n2=path[i+1]
                weight = self._g.get_edge_data(n1,n2)['weight']
                if weight==0:
                    weight=1
                edges_in_path.add((n1,n2,weight))
            # Add the new path and all edge weights
            for node in path:
                g.add_weighted_edges_from(edges_in_path)
        return g
    
    def _subtree(self, root):
        # Find paths to targets to build a new tree
        g = self._subgraph(root)
        #for node in g:
        #    print self.tree.node(node).data.taxon

        sub = Tree(weight=self.tree.weight, 
                   rooted=self.tree.rooted, 
                   name=self.tree.name, 
                   data=self.tree.dataclass, 
                   max_support=self.tree.max_support)

        sub.node(sub.root).data = self.tree.node(root).data
        def _add(old_id,new_id):
            """
            Walk from this node, using the id from the old tree, and the id
            from the new tree to both load the data from the old tree and link 
            to the correct node in the new tree.
            """
            for old_succ in g.successors_iter(old_id):
                to_add = Node(data = self.tree.node(old_succ).data)
                new_succ = sub.add(to_add,new_id)
                _add(old_succ,new_succ)
        _add(root,sub.root)

        # Delete nodes from old tree
        for node in g:
            #print "collapsing node",node
            self.collapse(self.tree,node)
            
        return sub
            
    def collapse(self, tree, id):
        """Deletes node from chain and relinks successors to predecessor: collapse(self, id)."""
        if id not in tree.chain:
            raise ChainException('Unknown ID: '+str(id))
        
        if not id==tree.root:
            prev_id=tree.chain[id].get_prev()
            tree.chain[prev_id].remove_succ(id)
            succ_ids=tree.chain[id].get_succ()
            for i in succ_ids:
                tree.chain[i].set_prev(prev_id)
            tree.chain[prev_id].add_succ(succ_ids)
            node=tree.chain[id]
            tree.kill(id)
            return node
        # The root needs to be weirdly joined up.
        # Pull a successor up and make it the root
        else:
            node = tree.chain[id]
            succ = tree.chain[id].get_succ()
            if len(succ)>0:
                tree.root = succ[0]
                tree.chain[tree.root].set_prev(None)
                if len(succ)>1:
                    for other_succ in succ[1:]:
                        tree.chain[succ[0]].add_succ(other_succ)
                        tree.chain[other_succ].set_prev(succ[0])
            else:
                tree.root = None
            tree.kill(id)
            #print "We won't collapse the root"
            return node
        
        
class TreeSplitter(object):

    def __init__(self, tree, max_size=40, annotater=UnrootedShortestPath):
        self.tree = tree
        self.max_size = max_size
        self.annotater = annotater(self.tree, max_size=max_size)
        self.annotater.annotate()
    
    def subtrees(self):
        """
        @return: generator of subtrees meeting size requirements.
        """
        #print "max_size",max_size
        sorted_nodes = self.annotater.sorted()
        # Keep track of the number removed to subtract from node totals after removing a subtree
        count_terminals = lambda node: self.tree.count_terminals(node)
        queue = sorted_nodes
        while len(queue)>0:
            #print "Queue",len(queue)
            prev = None
            cost = sorted(zip(map(self.annotater.cost,queue),queue),cmp=lambda x,y: cmp(x[0],y[0]))
            print [(self.tree.node(node).data.taxon,c) for c,node in cost]
            lowest_cost, best_root = cost[0]
            subtree = self.annotater.subtree(best_root)
            yield subtree
            self.annotater.annotate()
            queue = self.annotater.sorted()
        #yield self.tree
                


def unlink(tree, node):
    successors = list(tree.node(node).succ)
    #print "successors",successors
    for succ in successors:
        unlink(tree,succ)
    #print "unlinking",node
    tree.unlink(node)
    tree.kill(node)

if __name__=="__main__":
    from Bio.Nexus.Trees import Tree
    with open("/Users/patrick/experiment/amnh/subtree/just_mads.arab.branch") as handle:
        t = Tree(handle.read())
    n = 40
    s = TreeSplitter(t,max_size=n)
    
    taxa = set()
    terminals = 0
    s.tree.display()
    for tree in s.subtrees():
        tree.display()
        tax = tree.get_taxa()
        terminals+=len(tax)
        taxa.update(set(tax))
        #print tree
        print "Taxa",len(tax), tax
        print "terminals",len(tree.get_terminals())
        assert len(tax)<=n
    print "Terminals",terminals
    print "Taxa",len(taxa)
    assert len(taxa)==terminals
#
#def bifurcate(id):
#    node = tree.node(id)
#    succ = node.get_succ()
#    if not len(succ)>2:
#        return
#    keep = random.sample(succ,2)
#    split = keep[0]
#    for down_id in [i for i in succ if not i in keep]:
#        down_node = tree.node(down_id)
#        down_node.set_prev(split)
#        node.remove_succ(down_id)
#        tree.node(split).add_succ(down_id)
#    

