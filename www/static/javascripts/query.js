import { graph } from "./graph.js";

let searchButton = document.getElementById("search-button");
let searchInput = document.getElementById("search-input");
let searchDepth = document.getElementById("search-depth");

function searchCitationTree(pmid, depth) {
    const url = "/api/map/" + pmid;
    fetch(url, {
        depth: depth
    })
    .then(response => response.json())
    .then(data => {
            graph(data);
        });
}

function launchSearch() {
    let pmid = searchInput.value;
    let depth = searchDepth.value;
    searchCitationTree(pmid, depth);
}

searchButton.addEventListener("click", launchSearch)
searchInput.addEventListener("keyup", function(event) {
    // Enter key, not deprecated
    if (event.key === "Enter") {
        event.preventDefault();
        launchSearch();
    }
});


