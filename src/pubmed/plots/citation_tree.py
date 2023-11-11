import networkx as nx

def citation_tree_layer(G):
    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
    nx.draw(G, pos, with_labels=False, node_size=2, edge_color="black", alpha=0.75, width=0.2)
