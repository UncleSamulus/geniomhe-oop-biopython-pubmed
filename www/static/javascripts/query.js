import { graph } from "./graph.js"

const url = "/api/map/30917603";

fetch(url)
    .then(response => response.json())
    .then(data => {
        graph(data);
    })