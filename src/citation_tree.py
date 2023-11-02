#!/usr/bin/env python3
"""
Retrieve citation tree from a given paper
"""

from pubmed import explorer
import networkx as nx
import queue
from Bio import Entrez

Entrez.email = "samuel.ortion@etud.univ-evry.fr"

pmid = "37615633"
max_depth = 25

explorer = explorer.Explorer()

root = explorer.breadth_first_exploration(pmid, max_depth=max_depth, max_iterations=100)

G = nx.DiGraph()
q = queue.Queue()
current = root
q.put(current)
while not q.empty():
    current = q.get()
    for child in current.children:
        q.put(child)
        G.add_edge(current.data["pmid"], child.data["pmid"])

nx.write_graphml(G, f"tmp/citation_tree_{pmid}_d{max_depth}.graphml")
