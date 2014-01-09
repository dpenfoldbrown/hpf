'''
Created on Aug 18, 2009
@author: Patrick
'''
import networkx as nx
from hpf.function.term import Terms, Term
DEFAULT_QUERY="select t.id,t.acc from term t where t.acc in (%s)"

   
def _pydot(node):
    return Term(node.acc.replace("GO:",""),node.name)

class GODAG(object):
    '''
    '''

    def __init__(self, nodes=None,ebunch=None):
        self._nodes = []
        self._edges = []
        if nodes:
            self.add_nodes(nodes)
        if ebunch:
            self.add_edges(ebunch)

    def add_node(self, *nodes):
       for node in nodes:
            if hasattr(node,'__iter__'):
                self.add_nodes(node)
            else:
                self._nodes.append(node)
        
    def add_nodes(self, nodes):
        """
        @param nbunch: A container of Terms. 
        """
        for node in nodes:
            self.add_node(node)

    def add_edges(self, ebunch):
        """
        @param ebunch: A container of (Term1,Term2) tuples. 
        """
        self._edges += ebunch
        
    def load_edges(self,connection,
                   query=DEFAULT_QUERY):
        if not len(self._nodes) > 0:
            return
        terms_str = ",".join(["'%s'"%t.get_id() for t in self._nodes])
        query = query % (terms_str)
        cursor = connection.cursor()
        cursor.execute("drop table if exists terms1")
        cursor.execute("drop table if exists terms2")
        q = "create temporary table terms1 "+query
        #print q
        cursor.execute(q)
        q = "create temporary table terms2 "+query
        #print q
        cursor.execute(q)
        q = """select t1.acc,t2.acc from terms1 t1 join graph_path p
            on t1.id=p.term1_id and p.distance=1
            join terms2 t2 on p.term2_id=t2.id
        """
        #print q
        cursor.execute(q)
        edges = [(Term(t1),Term(t2)) for t1,t2 in cursor.fetchall()]
        cursor.close()
        self.add_edges(edges)

    def nodes(self):
        return self._nodes

    def graph(self):
        graph = nx.DiGraph()
        for node in self._nodes:
            graph.add_node(_pydot(node))
        for t1,t2 in self._edges:
            graph.add_edge(_pydot(t1),_pydot(t2))
        return graph
