import networkx as nx
import matplotlib.pyplot as plt
import click

@click.command()
@click.argument("input", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
def citation_tree(input, output):
    G = nx.read_graphml(input)
    citation_tree_layer(G)
    plt.savefig(output, dpi=600)

def citation_tree_layer(G):
    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
    nx.draw(G, pos, with_labels=False, node_size=2, edge_color="black", alpha=0.75, width=0.2)

if __name__ == "__main__":
    citation_tree()