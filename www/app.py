import json
from dotenv import load_dotenv, dotenv_values

from flask import Flask
from flask_caching import Cache
from flask import render_template
from flask import request

import networkx as nx
import pandas as pd
from networkx.readwrite import json_graph

from pubmed.explorer import Explorer

from Bio import Entrez

Entrez.email = "samuel.ortion@etud.univ-evry.fr"
explorer = Explorer()

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
    # Get depth from query string or default to 1
    depth = int(request.args.get("depth", 1))
    app.logger.info("Breadth first searching PMID:%s", pmid)
    root = explorer.breadth_first_exploration(pmid, max_depth=depth, max_iterations=1_000)
    app.logger.info("Done breadth first search")
    G = explorer.root_to_digraph(root)
    data = json_graph.node_link_data(G)
    return data


@app.route("/api/map/test", methods=["GET"])
def get_citation_map_test():
    with open("static/data/37847638.json", "r") as f:
        data = json.load(f)
    return data

@app.route("/api/trend/test", methods=["GET"])
def get_trend_test():
    data = pd.read_csv("../tmp/pubmed_cancer_2010-2022.csv")
    return data["date"].to_json()

@app.route("/api/search", methods=["GET"])
def query_keyword():
    keyword = request.args.get("keyword")
    start_date = "2010"
    end_date = "2022"
    data = explorer.search_info_csv(keyword, start_date, end_date)
    return data


if __name__=='__main__':
    app.run(debug=True)