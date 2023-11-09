import { citationGraph, trendGraph } from "./graph.js";

let searchButton = document.getElementById("search-button");
let searchInput = document.getElementById("search-input");
let searchDepth = document.getElementById("search-depth");


function searchKeyword(keyword) {
    const url = "/api/search/";
    fetch(url, {
        q: keyword
    
    })
    .then(
        response => response.json()
    )
    .then(data => {
            console.log(data);
    });
}

function searchCitationTree(pmid, depth) {
    const url = "/api/map/" + pmid;
    fetch(url, {
        depth: depth
    })
    .then(response => response.json())
    .then(data => {
            citationGraph(data);
        });
}


function searchTrend(pmid) {
    const url = "/api/trend/test";
    fetch(url)
    .then(response => response.json())
    .then(data => {
            trendGraph(data);
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


// (function() {
//     searchTrend();
// })();