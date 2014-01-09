from hpf.function.term import Term
from hpf.function.metric import Probability
from hpf.go import GODAG
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import networkx as nx


def plotIt():
    prob = Probability()
    prob.set_metric("a.1.1","A",1.0)
    prob.set_metric("a.1.1","B",0.83)
    prob.set_metric("a.1.1","C",0.84)
    prob.set_metric("a.1.1","D",0.1)
    prob.set_metric("a.1.1","E",0.01)
    prob.set_metric("a.1.1","Z",0.01)
    prob.set_metric("a.1.1","Y",0.01)
    prob.set_metric("a.1.1","X",0.01)
    prob.set_metric("a.1.1","W",0.01)
    prob.set_metric("a.1.1","V",0.01)
    prob.set_metric("a.1.1","U",0.01)
    prob.set_metric("a.1.1","T",0.01)
    prob.set_metric("a.1.1","S",0.01)
    prob.set_metric("a.1.1","R",0.01)
    prob.set_metric("a.1.1","Q",0.01)
    prob.set_metric("a.1.1","P",0.01)
    prob.set_metric("A",None,1.0)
    prob.set_metric("B",None,0.8)
    prob.set_metric("C",None,0.8)
    prob.set_metric("D",None,0.1)
    prob.set_metric("E",None,0.01)
    prob.set_metric("Z",None,1.0)
    prob.set_metric("Y",None,0.9)
    prob.set_metric("X",None,0.8)
    prob.set_metric("W",None,0.7)
    prob.set_metric("V",None,0.6)
    prob.set_metric("U",None,0.5)
    prob.set_metric("T",None,0.4)
    prob.set_metric("S",None,0.3)
    prob.set_metric("R",None,0.2)
    prob.set_metric("Q",None,0.1)
    prob.set_metric("P",None,0.0)

    nbunch = [(Term(m,"name"),prob.get_metric(sf,m)) for sf,m in prob if sf != None and m!= None]
    ebunch = [(Term("A"),Term("B")),
              (Term("A"),Term("C")),
              (Term("C"),Term("D")),
              (Term("D"),Term("E")),
              (Term("Z"),Term("Y")),
              (Term("Y"),Term("X")),
              (Term("X"),Term("W")),
              (Term("W"),Term("V")),
              (Term("V"),Term("U")),
              (Term("U"),Term("T")),
              (Term("T"),Term("S")),
              (Term("S"),Term("R")),
              (Term("R"),Term("Q")),
              (Term("Q"),Term("P")),]
    dag = GODAG(nbunch=nbunch,ebunch=ebunch)
    print dag.graph.nodes()
    matplotlib.pyplot.figure(figsize=(5, 5))
    # with nodes colored by degree sized by population
    node_color=[float(prob.get_metric(node.get_id(),None)) for node in dag.graph.nodes()]
    pos = nx.pydot_layout(dag.graph,prog="dot")
    labels = {}
    for node in dag.graph.nodes():
        labels[node] = node.acc+"\n"+node.name
    nx.draw(dag.graph,
            font_size=8,
            labels=labels,
            pos=pos,
            cmap="autumn",
            node_color=node_color,
            node_size=[50+prob.get_metric("a.1.1",node)*1000 for node, data in dag.graph.nodes(data=True)],
            #node_size=[100 for node, data in dag.graph.nodes(data=True)],
            with_labels=True)
    # scale the axes equally
    #matplotlib.pyplot.xlim(-5000, 500)
    #matplotlib.pyplot.ylim(-2000, 3500)
    filename = "/home/patrick/Sites/testGO.png"
    matplotlib.pyplot.savefig(filename)
