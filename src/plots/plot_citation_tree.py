import networkx as nx
import matplotlib.pyplot as plt

G = nx.read_graphml("tmp/citation_tree_37615633_d25.graphml")

pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
nx.draw(G, pos, with_labels=False, node_size=2, edge_color="black", alpha=0.75, width=0.2)

plt.savefig("./tmp/citation_tree_37615633_d25.png", dpi=600)
