import json
from dotenv import load_dotenv, dotenv_values

from flask import Flask
from flask_caching import Cache
from flask import render_template

import networkx as nx
from networkx.readwrite import json_graph

from pubmed.explorer import Explorer
explorer = Explorer()

from Bio import Entrez
Entrez.email = "samuel.ortion@etud.univ-evry.fr"

load_dotenv()
CONFIG = { **dotenv_values(".env") }
app = Flask(__name__)
app.config.from_mapping(CONFIG)
cache = Cache(app)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/test")
def test():
    return render_template("test.html")


@app.route("/api/map/<pmid>", methods=["GET"])
@cache.cached(timeout=50)
def get_citation_map(pmid):
    DEPTH = 3
    app.logger.info("Breadth first searching PMID:%s", pmid)
    root = explorer.breadth_first_exploration(pmid, max_depth=DEPTH, max_iterations=1_000)
    app.logger.info("Done breadth first search")
    G = explorer.root_to_digraph(root)
    data = json_graph.node_link_data(G)
    return data


@app.route("/api/map/test", methods=["GET"])
def get_citation_map_test():
    with open("static/data/37847638.json", "r") as f:
        data = json.load(f)
    return data


if __name__=='__main__':
    app.run(debug=True)