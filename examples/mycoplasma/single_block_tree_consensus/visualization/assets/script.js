//consensuses_colors = ['#ff8000', '#602870', '#983352','yellow', 'purple', 'orange', 'brown']

$(document).ready(function() {
    function show_info(info){
        name_p = document.getElementById('name');
        name_p.innerHTML += info.name;

        running_time_p = document.getElementById('running_time');
        running_time_p.innerHTML += info.running_time;

        sources_count_p = document.getElementById('sources_count');
        sources_count_p.innerHTML += info.sources_count;

        nodes_count_p = document.getElementById('poagraphs_count');
        nodes_count_p.innerHTML += info.poagraphs.length;
    }

    function show_blocks(data){
        var colors = d3.scaleOrdinal(d3.schemeCategory10);

        var svg = d3.selectAll("svg").filter("#blocks_svg"),
            width = +svg.attr("width"),
            height = +svg.attr("height"),
            node,
            link;

        svg.append('defs').append('marker')
            .attrs({'id':'arrowhead',
                'viewBox':'-0 -5 10 10',
                'refX':13,
                'refY':0,
                'orient':'auto',
                'markerWidth':13,
                'markerHeight':13,
                'xoverflow':'visible'})
            .append('svg:path')
            .attr('d', 'M 0,-5 L 10 ,0 L 0,5')
            .attr('fill', '#999')
            .style('stroke','none');

        var simulation = d3.forceSimulation()
            .force("link", d3.forceLink().id(function (d) {return d.id;}).distance(300)   .strength(1))
            .force("charge", d3.forceManyBody())
            .force("center", d3.forceCenter(width / 2, height / 2));

        d3.json("blocks.json", function (error, graph) {
            if (error) throw error;
            update(graph.edges, graph.nodes);
        })

        function update(links, nodes) {
            link = svg.selectAll(".link")
                .data(links)
                .enter()
                .append("line")
                .attr("class", "link")
                .attr('marker-end','url(#arrowhead)')
                .attr("stroke-dasharray", d=> d.active ? "0" : "10,10")

            link.append("title")
                .text(function (d) {return d.weight;});

            edgepaths = svg.selectAll(".edgepath")
                .data(links)
                .enter()
                .append('path')
                .attrs({
                    'class': 'edgepath',
                    'fill-opacity': 0,
                    'stroke-opacity': 0,
                    'id': function (d, i) {return 'edgepath' + i}
                })
                .style("pointer-events", "none");

            edgelabels = svg.selectAll(".edgelabel")
                .data(links)
                .enter()
                .append('text')
                .style("pointer-events", "none")
                .attrs({
                    'class': 'edgelabel',
                    'id': function (d, i) {return 'edgelabel' + i},
                    'font-size': 20,
                    'fill': '#333333'
                });

            edgelabels.append('textPath')
                .attr('xlink:href', function (d, i) {return '#edgepath' + i})
                .style("text-anchor", "middle")
                .style("pointer-events", "none")
                .attr("startOffset", "50%")
                .text(function (d) {return d.weight});

            node = svg.selectAll(".node")
                .data(nodes)
                .enter()
                .append("g")
                .attr("class", "node")
                .call(d3.drag()
                        .on("start", dragstarted)
                        .on("drag", dragged)
                        //.on("end", dragended)
                );

            node.append("circle")
                .attr("r", 25)
                .style("fill", function (d, i) {return colors(i);})

            node.append("title")
                .text(function (d) {return d.id;});

            node.append("text")
                .attr("dy", -3)
                .text(function (d) {return d.id;});

            simulation
                .nodes(nodes)
                .on("tick", ticked);

            simulation.force("link")
                .links(links);
        }

        function ticked() {
            link
                .attr("x1", function (d) {return d.source.x;})
                .attr("y1", function (d) {return d.source.y;})
                .attr("x2", function (d) {return d.target.x;})
                .attr("y2", function (d) {return d.target.y;});

            node
                .attr("transform", function (d) {return "translate(" + d.x + ", " + d.y + ")";});

            edgepaths.attr('d', function (d) {
                return 'M ' + d.source.x + ' ' + d.source.y + ' L ' + d.target.x + ' ' + d.target.y;
            });

            edgelabels.attr('transform', function (d) {
                if (d.target.x < d.source.x) {
                    var bbox = this.getBBox();

                    rx = bbox.x + bbox.width / 2;
                    ry = bbox.y + bbox.height / 2;
                    return 'rotate(180 ' + rx + ' ' + ry + ')';
                }
                else {
                    return 'rotate(0)';
                }
            });
        }

        function dragstarted(d) {
            if (!d3.event.active) simulation.alphaTarget(0.3).restart()
            d.fx = d.x;
            d.fy = d.y;
        }

        function dragged(d) {
            d.fx = d3.event.x;
            d.fy = d3.event.y;
        }

        function dragended(d) {
            if (!d3.event.active) simulation.alphaTarget(0);
            d.fx = undefined;
            d.fy = undefined;
        }
    }

    function show_sources_tables(poagraph_name, level_index){
        function draw_table(consensuses, sources){
            console.log(consensuses, sources);
            var sources_info_div = document.getElementById("sources_table");
            sources_info_div.innerHTML = "";
            var consensus_info_div = document.createElement('div')

            consensus_info_div.setAttribute("class", "consensus_info")

            var table = document.createElement('table');
            table.className = "sortable";
            var header  = table.createTHead();
            var header_row = header.insertRow(0);
            var head_id = header_row.insertCell(0);
            var head_name = header_row.insertCell(1);
            var head_title = header_row.insertCell(2);
            var head_group = header_row.insertCell(3);
            var head_bundle_id = header_row.insertCell(4);
            var consensuses_count = consensuses.length;
            head_id.innerHTML = "ID"
            head_name.innerHTML = "Name"
            head_title.innerHTML = "Title"
            head_group.innerHTML = "Group"
            head_bundle_id.innerHTML = "Bundle ID"

            for(var i=0; i< consensuses.length; i++){
                var cons_compatibility = header_row.insertCell(5+i);
                cons_compatibility.innerHTML = consensuses[i].ID;
                //cons_compatibility.innerHTML = consensuses[i].name;

            }
            body = document.createElement('tbody');

            for(var j=0;j<sources.length;j++){
                var row = body.insertRow(j);
                var id = row.insertCell(0);
                var name = row.insertCell(1);
                var title = row.insertCell(2);
                var group =row.insertCell(3)
                var bundle_id = row.insertCell(4);
                id.innerHTML = j;
                name.innerHTML = sources[j]['name'];
                title.innerHTML = sources[j]['title'];
                group.innerHTML = (sources[j]['group'].length) ? sources[j]['group'] : 'not defined';

                var sorted_levels_keys = Object.keys(sources[j]['bundle_ID']).sort()
                for(var k=0; k<sorted_levels_keys.length ; k++){
                    if(sorted_levels_keys[k] > level_index){
                        bundle_id.innerHTML = sources[j]['bundle_ID'][sorted_levels_keys[k]];
                        break;
                    }
                }


                for(var k=0;k<consensuses.length; k++){
                    var compatibility = row.insertCell(5+k);
                    c = consensuses[k]['sources_compatibility'][j];
                    if(consensuses[k].sources.indexOf(j) != -1){
                        compatibility.bgColor = "yellow"
                    }
                    compatibility.innerHTML = c;
                }
            }
            table.appendChild(body)
            table.createTFoot()
            consensus_info_div.appendChild(table)
            sources_info_div.appendChild(consensus_info_div);
            sorttable.makeSortable(table);

        }

        function filter_consensuses(consensuses_data){
            var consensuses_for_table = []
            var closed = []
            for(var i=0;i<consensuses_data.length; i++){
                var c = consensuses_data[i];
                if(c.level > level_index)
                {
                    var cc = 0;
                    for(cc=0; cc<closed.length; cc++)
                    {
                        if(closed[cc].ID == c.parent)
                        {
                            closed.push(c);
                            break;
                        }
                    }
                    if(cc == closed.length)
                    {
                        consensuses_for_table.push(c);
                        closed.push(c);
                    }
                }
                else
                {

                }
            }
            return consensuses_for_table;
        }

        consensuses_json_path = poagraph_name + "/consensuses.json"
        sources_json_path = poagraph_name + "/sources.json"
        $.getJSON(consensuses_json_path, function(consensuses_data) {
            consensuses_to_use = filter_consensuses(consensuses_data)
            $.getJSON(sources_json_path, function(sources_data){
                draw_table(consensuses_to_use, sources_data);
            });
        });


    }

    function setup_slider(poagraph_name){
        var slider = document.getElementById("slider");
        var output = document.getElementById("slider_value");
        output.innerHTML = slider.value;
        slider.oninput = function() {
            output.innerHTML = this.value;
        }

        var show_table_btn = document.getElementById("show_table_btn");
        show_table_btn.onclick = function() {
            show_sources_tables(poagraph_name, slider.value);
        }
    }

    function draw(consensuses){
        function transorm_to_tree_data(consensuses){
            function get_consensus_node(id){
                var c = consensuses[id];
                if(c.children.length == 0){
                    return {"name": "ID" + String(c.ID) + " : " + String(c.level) + " : " + String(c.sources),
                            "parent":c.parent,
                            "children":[]}
                }
                var children = []
                for(var j=0; j<c.children.length; j++){
                    children.push(get_consensus_node(c.children[j]))
                }
                return {"name": "ID" + String(c.ID) + " : " + String(c.level) + " : " + String(c.sources.length) + " sources",
                        "parent":c.parent,
                        "children":children}
            }

            treeData = []
            for(var i=1; i< consensuses.length; i++){
                if(consensuses[i].parent == 0){
                    treeData.push(get_consensus_node(i))
                }
            }
            return {"name": "",
                    "parent": "null",
                    "children": treeData};
        }

        var treeData = transorm_to_tree_data(consensuses);
                // Set the dimensions and margins of the diagram
        var margin = {top: 20, right: 90, bottom: 30, left: 90},
            width = 1960 - margin.left - margin.right,
            height = 500 - margin.top - margin.bottom;

        // append the svg object to the body of the page
        // appends a 'group' element to 'svg'
        // moves the 'group' element to the top left margin
        var svg = d3.select("#consensuses_tree").append("svg")
            .attr("width", width + margin.right + margin.left)
            .attr("height", height + margin.top + margin.bottom)
          .append("g")
            .attr("transform", "translate("
                  + margin.left + "," + margin.top + ")");

        var i = 0,
            duration = 750,
            root;

        // declares a tree layout and assigns the size
        var treemap = d3.tree().size([height, width]);

        // Assigns parent, children, height, depth
        root = d3.hierarchy(treeData, function(d) { return d.children; });
        root.x0 = height / 2;
        root.y0 = 0;

        // Collapse after the second level
        root.children.forEach(collapse);

        update(root);

        // Collapse the node and all it's children
        function collapse(d) {
          if(d.children) {
            d._children = d.children
            d._children.forEach(collapse)
            d.children = null
          }
        }

        function update(source) {

          // Assigns the x and y position for the nodes
          var treeData = treemap(root);

          // Compute the new tree layout.
          var nodes = treeData.descendants(),
              links = treeData.descendants().slice(1);

          // Normalize for fixed-depth.
          nodes.forEach(function(d){ d.y = d.depth * 180});

          // ****************** Nodes section ***************************

          // Update the nodes...
          var node = svg.selectAll('g.node')
              .data(nodes, function(d) {return d.id || (d.id = ++i); });

          // Enter any new modes at the parent's previous position.
          var nodeEnter = node.enter().append('g')
              .attr('class', 'node')
              .attr("transform", function(d) {
                return "translate(" + source.y0 + "," + source.x0 + ")";
            })
            .on('click', click);

          // Add Circle for the nodes
          nodeEnter.append('circle')
              .attr('class', 'node')
              .attr('r', 1e-6)
              .style("fill", function(d) {
                  return d._children ? "lightsteelblue" : "#fff";
              });

          // Add labels for the nodes
          nodeEnter.append('text')
              .attr("dy", ".35em")
              .attr("x", function(d) {
                  return d.children || d._children ? -13 : 13;
              })
              .attr("text-anchor", function(d) {
                  return d.children || d._children ? "end" : "start";
              })
              .text(function(d) { return d.data.name; });

          // UPDATE
          var nodeUpdate = nodeEnter.merge(node);

          // Transition to the proper position for the node
          nodeUpdate.transition()
            .duration(duration)
            .attr("transform", function(d) {
                return "translate(" + d.y + "," + d.x + ")";
             });

          // Update the node attributes and style
          nodeUpdate.select('circle.node')
            .attr('r', 10)
            .style("fill", function(d) {
                return d._children ? "lightsteelblue" : "#fff";
            })
            .attr('cursor', 'pointer');


          // Remove any exiting nodes
          var nodeExit = node.exit().transition()
              .duration(duration)
              .attr("transform", function(d) {
                  return "translate(" + source.y + "," + source.x + ")";
              })
              .remove();

          // On exit reduce the node circles size to 0
          nodeExit.select('circle')
            .attr('r', 1e-6);

          // On exit reduce the opacity of text labels
          nodeExit.select('text')
            .style('fill-opacity', 1e-6);

          // ****************** links section ***************************

          // Update the links...
          var link = svg.selectAll('path.link')
              .data(links, function(d) { return d.id; });

          // Enter any new links at the parent's previous position.
          var linkEnter = link.enter().insert('path', "g")
              .attr("class", "link")
              .attr('d', function(d){
                var o = {x: source.x0, y: source.y0}
                return diagonal(o, o)
              });

          // UPDATE
          var linkUpdate = linkEnter.merge(link);

          // Transition back to the parent element position
          linkUpdate.transition()
              .duration(duration)
              .attr('d', function(d){ return diagonal(d, d.parent) });

          // Remove any exiting links
          var linkExit = link.exit().transition()
              .duration(duration)
              .attr('d', function(d) {
                var o = {x: source.x, y: source.y}
                return diagonal(o, o)
              })
              .remove();

          // Store the old positions for transition.
          nodes.forEach(function(d){
            d.x0 = d.x;
            d.y0 = d.y;
          });

          // Creates a curved (diagonal) path from parent to the child nodes
          function diagonal(s, d) {

            path = `M ${s.y} ${s.x}
                    C ${(s.y + d.y) / 2} ${s.x},
                      ${(s.y + d.y) / 2} ${d.x},
                      ${d.y} ${d.x}`

            return path
          }

          // Toggle children on click.
          function click(d) {
            if (d.children) {
                d._children = d.children;
                d.children = null;
              } else {
                d.children = d._children;
                d._children = null;
              }
            update(d);
          }
        }
    }

    $.getJSON("info.json", function(info) {
      show_info(info);
      var poagraphs_names = info.poagraphs;
      var sources_json_path = poagraphs_names[0] + "/sources.json";

      setup_slider(poagraphs_names[0]);

      var consensuses_json_path = poagraphs_names[0] + "/consensuses.json"
      $.getJSON(consensuses_json_path, function(consensuses) {
          draw(consensuses);
      });

    });

    $.getJSON("blocks.json", function(json) {
      show_blocks(json);
    });

    //add_consensuses_options();
    //draw_visualization();
});


function show_blocks_cytoscape(blocks){
    cytoscape({
        container: document.getElementById('sep_blocks_info'),
        elements: blocks,
        style: [
          {
            selector: 'node',
            style: {
              'label': 'data(id)',
              'width': '15px',
              'height': '15px',
              'font-size': '8px',
              'text-halign': 'center',
              'text-valign': 'center',
              'background-color' : 'white'
            }
          }, {
            selector: 'edge',
            style: {
              'curve-style': 'bezier',
              'label': get_label,
              'line-color': get_edge_label_color,
              'width': '3px',
              'target-arrow-shape': 'triangle',
              'target-arrow-color': get_edge_label_color
          }
        }
        ],
        layout: {
          name: 'breadthfirst'
        }
  });

  function get_label(ele) {
    srcIDs = ele.data('srcID')
    return String(srcIDs.length)
//    if(srcIDs.length > 5)
//        return '';
//    return srcIDs;
    }

  function get_edge_label_color(ele) {
    return 'Red';
    srcID = ele.data('srcID');
    cl = [  'Red',
            'Green'	,
            'Yellow'	,
            'Blue'	,
            'Orange'	,
            'Purple'	,
            'Cyan'	,
            'Magenta'	,
            'Lime'	,
            'Pink'	,
            'Teal'	,
            'Lavender',
            'Brown'	,
            'Beige'	,
            'Maroon'	,
            'Mint'	,
            'Olive'	,
            'Coral'	,
            'Navy'	,
            'Grey'	,
            'White'	,
            'Black']
    return cl[srcID%cl.length];
  }

//  function get_edge_width(ele) {
//    edge_width = Math.abs(Math.log(ele.data('weight') * 20))
//    return edge_width.toString() + 'px';
//  }

//  function get_node_width(ele) {
//    node_width = Math.abs(Math.log(ele.data('weight') * 1000))
//    return node_width.toString() + 'px';
//  }
}

function add_consensuses_options(){
//    for(i=0; i<info.levels.length;i++)
//    {
//        let level_value = info.levels[i];
//        let level_index = i;
//        var radioBtn_p = $('<p></p>')
//
//        var radioBtn = $('<input type="radio" name="level" value="' + String(level_value) + '">');
//        radioBtn.click( function(){
//            show_sources_info_table(level_index);
//        });
//        radioBtn.appendTo(radioBtn_p);
//
//        var radioBtn_label = $('<label>' + String(level_value) + '</label>')
//        radioBtn_label.appendTo(radioBtn_p);
//
//        radioBtn_p.appendTo('#consensuses_info');
//    }
}

function draw_visualization(){
    var cy = cytoscape({
    container: document.getElementById('cy'),
    elements: poagraph,
    style: [
      {
        selector: 'node',
        style: {
          'label': 'data(nucleobase)',
          'width': get_node_width,
          'height': get_node_width,
          'font-size': '4px',
          'text-halign': 'center',
          'text-valign': 'center',
          'background-color' : 'white'
        }
      }, {
        selector: 'edge',
        style: {
          'curve-style': 'haystack',
          'label': 'data(title)',
          'line-color': get_edge_label_color, //to differentiate consensuses
          'width': get_edge_width, //to show popular edges
      }
    },{
      selector: '.aligned',
      style: {
        'line-style': 'dotted',
      }
    },{
      selector: '.consensus',
      style: {
        'curve-style': 'bezier',
      }
    }
    ],
    layout: {
      name: 'preset'
    }
  });

  function get_edge_label_color(ele) {

    consensus_id = ele.data('consensus');
    if(consensus_id === -1){
      return 'white';
    }
    return consensuses_colors[consensus_id%consensuses_colors.length];
  }

  function get_edge_width(ele) {
    edge_width = Math.abs(Math.log(ele.data('weight') * 20))
    return edge_width.toString() + 'px';
  }

  function get_node_width(ele) {
    node_width = Math.abs(Math.log(ele.data('weight') * 1000))
    return node_width.toString() + 'px';
  }

}

function show_consensuses_tree2(poagraph_name){
    var margin = {top: 20, right: 120, bottom: 20, left: 120},
	width = 1260 - margin.right - margin.left,
	height = 700 - margin.top - margin.bottom;

    var i = 0,
        duration = 750,
        root;

    var tree = d3.layout.tree()
        .size([height, width]);

    var diagonal = d3.svg.diagonal()
        .projection(function(d) { return [d.y, d.x]; });

//    var tree_svg = document.getElementById("tree_display");
    var svg = d3.select("#tree_display").attr("width", width + margin.right + margin.left)
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    //root = treeData[0];
    root = treeData[0];
    root.x0 = height / 2;
    root.y0 = 0;

    update(root);

    d3.select(self.frameElement).style("height", "500px");

    function update(source) {

      // Compute the new tree layout.
      var nodes = tree.nodes(root).reverse(),
          links = tree.links(nodes);

      // Normalize for fixed-depth.
      nodes.forEach(function(d) { d.y = d.depth * 180; });

      // Update the nodes…
      var node = svg.selectAll("g.node")
          .data(nodes, function(d) { return d.id || (d.id = ++i); });

      // Enter any new nodes at the parent's previous position.
//      var nodeEnter = node.enter().append("g")
//          .attr("class", "node")
//          .attr("transform", function(d) { return "translate(" + source.y0 + "," + source.x0 + ")"; })
//          .on("click", click);
      var nodeEnter = node.enter().append("g")
          .attr("class", "node")
          .attr("transform", function(d) { return "translate(" + source.y0 + "," + source.x0 + ")"; })
          .on("click", function() {
            click(this);
        });


      nodeEnter.append("circle")
          .attr("r", 1e-6)
          .style("fill", function(d) { return d._children ? "lightsteelblue" : "#fff"; });

      nodeEnter.append("text")
          .attr("x", function(d) { return d.children || d._children ? -13 : 13; })
          .attr("dy", ".35em")
          .attr("text-anchor", function(d) { return d.children || d._children ? "end" : "start"; })
          .text(function(d) { return (d.name + ': ' + d.val + (d.sources ? " sources: " + d.sources : "")); })
          .style("fill-opacity", 1e-6);

      // Transition nodes to their new position.
      var nodeUpdate = node.transition()
          .duration(duration)
          .attr("transform", function(d) { return "translate(" + d.y + "," + d.x + ")"; });

      nodeUpdate.select("circle")
          .attr("r", 10)
          .style("fill", function(d) { return d._children ? "lightsteelblue" : "#fff"; });

      nodeUpdate.select("text")
          .style("fill-opacity", 1);

      // Transition exiting nodes to the parent's new position.
      var nodeExit = node.exit().transition()
          .duration(duration)
          .attr("transform", function(d) { return "translate(" + source.y + "," + source.x + ")"; })
          .remove();

      nodeExit.select("circle")
          .attr("r", 1e-6);

      nodeExit.select("text")
          .style("fill-opacity", 1e-6);

      // Update the links…
      var link = svg.selectAll("path.link")
          .data(links, function(d) { return d.target.id; });

      // Enter any new links at the parent's previous position.
      link.enter().insert("path", "g")
          .attr("class", "link")
          .attr("d", function(d) {
            var o = {x: source.x0, y: source.y0};
            return diagonal({source: o, target: o});
          });

      // Transition links to their new position.
      link.transition()
          .duration(duration)
          .attr("d", diagonal);

      // Transition exiting nodes to the parent's new position.
      link.exit().transition()
          .duration(duration)
          .attr("d", function(d) {
            var o = {x: source.x, y: source.y};
            return diagonal({source: o, target: o});
          })
          .remove();

      // Stash the old positions for transition.
      nodes.forEach(function(d) {
        d.x0 = d.x;
        d.y0 = d.y;
      });
    }

    function click(d) {
      show_node_info_table(d.textContent)

      if (d.children) {
        d._children = d.children;
        d.children = null;
      } else {
        d.children = d._children;
        d._children = null;
      }
      update(d);

    }
}

function show_sources_info_table(level_index){
    var sources_info_div = document.getElementById("sources_info");
    sources_info_div.innerHTML = "";
    var consensus_info_div = document.createElement('div')

    consensus_info_div.setAttribute("class", "consensus_info")

    var table = document.createElement('table');
    table.className = "sortable";
    var header  = table.createTHead();
    var header_row = header.insertRow(0);
    var head_id = header_row.insertCell(0);
    var head_name = header_row.insertCell(1);
    var head_title = header_row.insertCell(2);
    var head_group = header_row.insertCell(3);
    var head_bundle_id = header_row.insertCell(4);
    var cses = get_consensuses_for_table(level_index);
    var consensuses_count = cses.length;
    head_id.innerHTML = "ID"
    head_name.innerHTML = "Name"
    head_title.innerHTML = "Title"
    head_group.innerHTML = "Group"
    head_bundle_id.innerHTML = "Bundle ID"

    for(var i=0; i< cses.length; i++){
        var cons_compatibility = header_row.insertCell(5+i);
        cons_compatibility.innerHTML = cses[i].name;
    }
    body = document.createElement('tbody')



    for(var j=0;j<sources.data.length;j++){
        var row = body.insertRow(j);
        var id = row.insertCell(0);
        var name = row.insertCell(1);
        var title = row.insertCell(2);
        var group =row.insertCell(3)
        var bundle_id = row.insertCell(4);
        id.innerHTML = j;
        name.innerHTML = sources.data[j]['name'];
        title.innerHTML = sources.data[j]['title'];
        group.innerHTML = (sources.data[j]['group'].length) ? sources.data[j]['group'] : 'not defined';

        var sorted_levels_keys = Object.keys(sources.data[j]['bundle_ID']).sort()
        for(var k=0; k<sorted_levels_keys.length ; k++){
            if(sorted_levels_keys[k] > level_index){
                bundle_id.innerHTML = sources.data[j]['bundle_ID'][sorted_levels_keys[k]];
                break;
            }
        }


        for(var k=0;k<cses.length; k++){
            var compatibility = row.insertCell(5+k);
            c = cses[k]['sources_compatibility'][j];
            if(cses[k].sources.indexOf(j) != -1){
                compatibility.bgColor = "yellow"
            }
            compatibility.innerHTML = c;
        }
    }
    table.appendChild(body)
    table.createTFoot()
    consensus_info_div.appendChild(table)
    sources_info_div.appendChild(consensus_info_div);
    sorttable.makeSortable(table);
}

function show_node_info_table(node_name){
    node_children = []
    if(node_name == "All sequences"){
        node_id = -1

        for(var i=0;i<consensuses.c.length;i++)
        {
            if(consensuses.c[i].parent == -1)
            {
                node_children.push(consensuses.c[i])
            }
        }
    }
    else{
        node_id = parseInt(node_name.substr(7));
        node_children_ids = consensuses.c[node_id].children;
        for(var i=0; i<node_children_ids.length;i++)
        {
            node_children.push(consensuses.c[node_children_ids[i]]);
        }
    }

    var sources_info_div = document.getElementById("node_info");
    sources_info_div.innerHTML = "";
    var consensus_info_div = document.createElement('div')
    consensus_info_div.setAttribute("class", "consensus_info")

    var table = document.createElement('table');
    table.className = "sortable";
    var header  = table.createTHead();
    var header_row = header.insertRow(0);
    var head_id = header_row.insertCell(0);
    var head_name = header_row.insertCell(1);
    var head_title = header_row.insertCell(2);
    var head_group = header_row.insertCell(3);
    var head_bundle_id = header_row.insertCell(4);
    var cses = node_children;
    var consensuses_count = cses.length;
    head_id.innerHTML = "ID"
    head_name.innerHTML = "Name"
    head_title.innerHTML = "Title"
    head_group.innerHTML = "Group"
    head_bundle_id.innerHTML = "Bundle ID"

    for(var i=0; i< cses.length; i++){
        var cons_compatibility = header_row.insertCell(5+i);
        cons_compatibility.innerHTML = cses[i].name;
    }
    body = document.createElement('tbody')


    var valid_source_id = 0;
    for(var j=0;j<sources.data.length;j++){
        if(node_id != -1 && consensuses.c[node_id].sources.indexOf(j) == -1){
            continue;
        }

        var row = body.insertRow(valid_source_id);
        var id = row.insertCell(0);
        var name = row.insertCell(1);
        var title = row.insertCell(2);
        var group =row.insertCell(3)
        var bundle_id = row.insertCell(4);
        id.innerHTML = j;
        name.innerHTML = sources.data[j]['name'];
        title.innerHTML = sources.data[j]['title'];
        group.innerHTML = (sources.data[j]['group'].length) ? sources.data[j]['group'] : 'not defined';

        var sorted_levels_keys = Object.keys(sources.data[j]['bundle_ID']).sort()
        if(node_id == -1)
        {
            var level = 0;
        }
        else
        {
            var level = consensuses.c[node_id].level
        }
        for(var k=0; k<sorted_levels_keys.length ; k++){
            if(sorted_levels_keys[k] > level){
                bundle_id.innerHTML = sources.data[j]['bundle_ID'][sorted_levels_keys[k]];
                break;
            }
        }

//        bundle_id.innerHTML = sources.data[j]['bundle_ID'][consensuses.c[node_id].level];

        for(var k=0;k<cses.length; k++){
            var compatibility = row.insertCell(5+k);
            c = cses[k]['sources_compatibility'][j];
            compatibility.innerHTML = c;
        }
        valid_source_id += 1;

    }
    table.appendChild(body)
    table.createTFoot()
    consensus_info_div.appendChild(table)
    sources_info_div.appendChild(consensus_info_div);
    sorttable.makeSortable(table);
}

function get_consensuses_for_table(level_index){
    var consensuses_for_table = []
    var closed = []
    for(var i=0;i<consensuses.c.length; i++){
        c = consensuses.c[i];
        if(c.level > level_index)
        {
            var cc = 0;
            for(cc=0; cc<closed.length; cc++)
            {
                if(closed[cc].id == c.parent)
                {
                    closed.push(c);
                    break;
                }
            }
            if(cc == closed.length)
            {
                consensuses_for_table.push(c);
                closed.push(c);
            }
        }
    }
    return consensuses_for_table;
}

