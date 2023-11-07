import { graph } from "./graph.js";

let searchButton = document.getElementById("search-button");
let searchInput = document.getElementById("search-input");
let searchDepth = document.getElementById("search-depth");

searchButton.addEventListener("click", () => {
    let pmid = searchInput.value;
    let depth = searchDepth.value;
    const url = "/api/map/" + pmid;
    fetch(url, {
        depth: depth
    })
    .then(response => response.json())
    .then(data => {
            graph(data);
        });
    });


