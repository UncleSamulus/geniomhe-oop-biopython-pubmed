import { citationGraph, trendGraph } from "./graph.js";

const pmid = "test";
const url = "/api/map/" + pmid;
    fetch(url, {
        depth: 2
    })
    .then(response => response.json())
    .then(data => {
            citationGraph(data);
        });
