function citationGraph(data) {

    let margin = {top: 10, right: 30, bottom: 30, left: 40},
        width = 100 - margin.left - margin.right,
        height = 100 - margin.top - margin.bottom;

    let svg = d3.select("#citation-graph")
        .append("svg")
        .append("g")
            .attr("transform",
                    "translate(" + margin.left + "," + margin.top + ")");


    let link = svg
        .selectAll("line")
        .data(data.links)
        .enter()
        .append("line")
        .style("stroke", "#aaa")
    
    let node = svg
        .selectAll("circle")
        .data(data.nodes)
        .enter()
        .append("circle")
            .attr("r", 10)
            .style("fill", "#69b3a2")

    let simulation = d3.forceSimulation(data.nodes)
        .force("link", d3.forceLink()
                .id(function(d) { return d.id; })
                .links(data.links)
            )
            .force("charge", d3.forceManyBody().strength(-400))
            .force("center", d3.forceCenter(width/2, height/2))
            .on("end", ticked);

    function ticked() {
        link
            .attr("x1", function(d) { return d.source.x; })
            .attr("y1", function(d) { return d.source.y; })
            .attr("x2", function(d) { return d.target.x; })
            .attr("y2", function(d) { return d.target.y; });
    
        node
            .attr("cx", function (d) { return d.x; })
            .attr("cy", function(d) { return d.y; });
        }
}

function trendGraph(timeseries) {
    let margin = {top: 10, right: 30, bottom: 30, left: 40},
        width = 100 - margin.left - margin.right,
        height = 100 - margin.top - margin.bottom;

    let svg = d3.select("#dataviz")
        .append("svg")
        .append("g")
            .attr("transform",
                    "translate(" + margin.left + "," + margin.top + ")");

    let x = d3.scaleTime()
        .domain(d3.extent(timeseries, function(d) { return d.date; }))
        .range([ 0, width ]);

    svg.append("g")
        .attr("transform", "translate(0," + height + ")")
        .call(d3.axisBottom(x));

    
    let y = d3.scaleLinear()
        .domain([0, d3.max(timeseries, function(d) { return +d.count; })])
        .range([ height, 0 ]);

    svg.append("g")
        .call(d3.axisLeft(y));

    svg.append("path")
        .datum(timeseries)
        .attr("fill", "none")
        .attr("stroke", "steelblue")
        .attr("stroke-width", 1.5)
        .attr("d", d3.line()
                .x(function(d) { return x(d.date) })
                .y(function(d) { return y(d.count) })
            )
}

export {
    citationGraph,
    trendGraph
}