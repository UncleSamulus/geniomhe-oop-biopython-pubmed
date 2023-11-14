#!/usr/bin/env python3
"""
Retrieve citation tree from a given paper
"""

from pubmed import explorer
import networkx as nx
import queue
from Bio import Entrez

Entrez.email = "samuel.ortion@etud.univ-evry.fr"

pmid = "32316364"
max_depth = 2
sense = "cited"

exp = explorer.Explorer()

root = exp.bidirectional_breadth_first_exploration(pmid, max_depth=max_depth)
G = exp.biroot_to_digraph(root)
nx.write_graphml(G, f"tmp/citation_tree_dual{pmid}_d{max_depth}.graphml")