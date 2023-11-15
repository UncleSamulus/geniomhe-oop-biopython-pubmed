import networkx as nx
import matplotlib.pyplot as plt
import click

@click.command()
@click.argument("citing", type=click.Path(exists=True))
@click.argument("cited", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
def citation_tree(citing, cited, output):
    G1 = nx.read_graphml(citing)
    G2 = nx.read_graphml(cited)
    G2_rev = G2.reverse()
    G = nx.compose(G1, G2_rev)
    citation_tree_layer(G)
    plt.savefig(output, dpi=600)

def citation_tree_layer(G):
    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
    nx.draw(G, pos, with_labels=False, node_size=0.2, edge_color="black", alpha=0.5, width=0.1, arrowstyle='-|>')

if __name__ == "__main__":
    citation_tree()