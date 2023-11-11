import queue

import click
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from Bio import Entrez
Entrez.email = "samuel.ortion@etud.univ-evry.fr"

import pubmed.explorer
from pubmed.plots import citation_tree

explorer = pubmed.explorer.Explorer()

@click.command()
@click.option("--pmid", "-p", help="Paper ID")
@click.option("--max-depth", "-d", default=5, help="Max depth of the tree")
def citation_map(pmid, max_depth=25, max_iterations=10_000):
    root = explorer.breadth_first_exploration(pmid, max_depth=max_depth, max_iterations=max_iterations)
    G = explorer.root_to_digraph(root)
    nx.write_graphml(G, f"tmp/citation_tree_{pmid}_d{max_depth}.graphml")
    citation_tree.citation_tree_layer(G)
    plt.savefig(f"tmp/citation_tree_{pmid}_d{max_depth}.png")


def main():
    citation_map()
    # start_date = 2010
    # end_date = 2022
    # query = "cancer"
    # explorer = pubmed.explorer.Explorer()
    # data = list(
    #     map(
    #         lambda item: explorer.extract_info(item[1], item[0]),
    #         explorer.query(query, start_date, end_date),
    #     )
    # )
    # if None in data:
    #     data.remove(None)
    # data_series = dict(
    #     pmid=[],
    #     title=[],
    #     doi=[],
    #     main_author=[],
    #     date=[],
    #     keywords=[],
    # )

    # for row in data:
    #     for key in row:
    #         data_series[key].append(row[key])

    # df = pd.DataFrame(data_series)
    # df.to_csv("tmp/test.csv")


if __name__ == "__main__":
    main()
