function get_density_graph(i){
    var density_graph = document.createElementNS(d3.namespaces.svg, 'svg');
        density_graph.setAttribute("class", "density_graph")
        density_graph.setAttribute("width", 560)
        density_graph.setAttribute("height",  250)

    var width = +density_graph.getAttribute("width"),
        height = +density_graph.getAttribute("height"),
        margin = {top: 20, right: 30, bottom: 30, left: 40};

    var x = d3.scaleLinear()
                .domain([0, 100])
                .range([margin.left, width - margin.right]);

    var y = d3.scaleLinear()
                .domain([0, 1])
                .range([height - margin.bottom, margin.top]);
    var piesek = d3.select(density_graph);

              piesek.append("text")
                    .attr("x", (width / 2))
                    .attr("y", 20)
                    .attr("text-anchor", "middle")
                    .text("Sequences compatibility to the consensus");
              piesek.append("g")
                    .attr("class", "axis axis--x")
                    .attr("transform", "translate(0," + (height - margin.bottom) + ")")
                    .call(d3.axisBottom(x))
                    .append("text")
                    .attr("x", width - margin.right)
                    .attr("y", -6)
                    .attr("fill", "#000")
                    .attr("text-anchor", "end")
                    .attr("font-weight", "bold")
                    .text("Compatibility");

              piesek.append("g")
                    .attr("class", "axis axis--y")
                    .attr("transform", "translate(" + margin.left + ",0)")
                    .call(d3.axisLeft(y).ticks(null, "%"));

    var myData = consensuses.data[i].sources_compatibility.map(function(x){return x*100;});

    var data_length = myData.length,
        bins = d3.histogram().domain(x.domain()).thresholds(40)(myData),
        density = kernelDensityEstimator(kernelEpanechnikov(7), x.ticks(40))(myData);
            piesek.insert("g", "*")
                  .attr("fill", "#bbb")
                  .selectAll("rect")
                  .data(bins)
                  .enter().append("rect")
                  .attr("x", function(d)
                   {
                    return x(d.x0) + 1;
                     })
                  .attr("y", function(d) {
                   return y(d.length / data_length);
                    })
                  .attr("width", function(d) {
                   return x(d.x1) - x(d.x0) - 1;
                    })
                  .attr("height", function(d) {
                   return y(0) - y(d.length / data_length);
                    });

            piesek.append("path")
                  .datum(density)
                  .attr("fill", "none")
                  .attr("stroke", "#000")
                  .attr("stroke-width", 1.5)
                  .attr("stroke-linejoin", "round")
                  .attr("d",  d3.line()
                  .curve(d3.curveBasis)
                  .x(function(d) { return x(d[0]); })
                  .y(function(d) { return y(d[1]); }));

    return piesek.node();
};

function kernelDensityEstimator(kernel, X) {
  return function(V) {
    return X.map(function(x) {
      return [x, d3.mean(V, function(v) { return kernel(x - v); })];
    });
  };
}

function kernelEpanechnikov(k) {
  return function(v) {
    return Math.abs(v /= k) <= 1 ? 0.75 * (1 - v * v) / k : 0;
  };
}