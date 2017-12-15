consensuses_colors = ['#ff8000', '#602870', '#983352','yellow', 'purple', 'orange', 'brown']

document.addEventListener("DOMContentLoaded", function() {
    add_poagraph_info();
    add_consensuses_options();
    draw_visualization();
    show_consensuses_tree();
    //show_sources_info_table();
    setup_slider();
});

function setup_slider(){
    var slider = document.getElementById("slider");
    var output = document.getElementById("slider_value");
    output.innerHTML = slider.value;
    slider.oninput = function() {
        output.innerHTML = this.value;
    }

    var show_table_btn = document.getElementById("show_table_btn");
    show_table_btn.onclick = function() {
        show_sources_info_table(slider.value);
    }
}

function add_poagraph_info(){
    name_p = document.getElementById('name');
    name_p.innerHTML += info.name;

    consensus_algorithm_p = document.getElementById('consensus_algorithm');
    consensus_algorithm_p.innerHTML += info.consensus_algorithm;

    running_time_p = document.getElementById('running_time');
    running_time_p.innerHTML += info.running_time;

    sources_count_p = document.getElementById('sources_count');
    sources_count_p.innerHTML += info.sources_count;

    nodes_count_p = document.getElementById('nodes_count');
    nodes_count_p.innerHTML += info.nodes_count;

    sequences_per_node_p = document.getElementById('sequences_per_node');
    sequences_per_node_p.innerHTML += info.sequences_per_node;
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

function show_consensuses_tree(){


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
      var nodeEnter = node.enter().append("g")
          .attr("class", "node")
          .attr("transform", function(d) { return "translate(" + source.y0 + "," + source.x0 + ")"; })
          .on("click", click);

      nodeEnter.append("circle")
          .attr("r", 1e-6)
          .style("fill", function(d) { return d._children ? "lightsteelblue" : "#fff"; });

      nodeEnter.append("text")
          .attr("x", function(d) { return d.children || d._children ? -13 : 13; })
          .attr("dy", ".35em")
          .attr("text-anchor", function(d) { return d.children || d._children ? "end" : "start"; })
          .text(function(d) { return d.name; })
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
      childrennames = []
      if(d.children)
      {
        for(var i; i<d.children.length;i++)
      {
        childrennames.push(d.children[i].name)
      }
      show_node_info_table(childrennames)
      } else
      {
      for(var i; i<d._children.length;i++)
      {
        childrennames.push(d._children[i].name)
      }
      show_node_info_table(childrennames)
      }

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
            compatibility.innerHTML = c;
        }
    }
    table.appendChild(body)
    table.createTFoot()
    consensus_info_div.appendChild(table)
    sources_info_div.appendChild(consensus_info_div);
    sorttable.makeSortable(table);
}

function show_node_info_table(ids){
    var sources_info_div = document.getElementById("node_info");
    sources_info_div.innerHTML = ids;
    console.log(ids);
    var consensus_info_div = document.createElement('div')
    consensus_info_div.setAttribute("class", "consensus_info")

//    var table = document.createElement('table');
//    table.className = "sortable";
//    var header  = table.createTHead();
//    var header_row = header.insertRow(0);
//    var head_id = header_row.insertCell(0);
//    var head_name = header_row.insertCell(1);
//    var head_title = header_row.insertCell(2);
//    var head_group = header_row.insertCell(3);
//    var head_bundle_id = header_row.insertCell(4);
//    var cses = get_consensuses_for_table(level_index);
//    var consensuses_count = cses.length;
//    head_id.innerHTML = "ID"
//    head_name.innerHTML = "Name"
//    head_title.innerHTML = "Title"
//    head_group.innerHTML = "Group"
//    head_bundle_id.innerHTML = "Bundle ID"
//
//    for(var i=0; i< cses.length; i++){
//        var cons_compatibility = header_row.insertCell(5+i);
//        cons_compatibility.innerHTML = cses[i].name;
//    }
//    body = document.createElement('tbody')
//
//
//
//    for(var j=0;j<sources.data.length;j++){
//        var row = body.insertRow(j);
//        var id = row.insertCell(0);
//        var name = row.insertCell(1);
//        var title = row.insertCell(2);
//        var group =row.insertCell(3)
//        var bundle_id = row.insertCell(4);
//        id.innerHTML = j;
//        name.innerHTML = sources.data[j]['name'];
//        title.innerHTML = sources.data[j]['title'];
//        group.innerHTML = (sources.data[j]['group'].length) ? sources.data[j]['group'] : 'not defined';
//
//        var sorted_levels_keys = Object.keys(sources.data[j]['bundle_ID']).sort()
//        for(var k=0; k<sorted_levels_keys.length ; k++){
//            if(sorted_levels_keys[k] > level_index){
//                bundle_id.innerHTML = sources.data[j]['bundle_ID'][sorted_levels_keys[k]];
//                break;
//            }
//        }
//
//
//        for(var k=0;k<cses.length; k++){
//            var compatibility = row.insertCell(5+k);
//            c = cses[k]['sources_compatibility'][j];
//            compatibility.innerHTML = c;
//        }
//    }
//    table.appendChild(body)
//    table.createTFoot()
//    consensus_info_div.appendChild(table)
    sources_info_div.appendChild(consensus_info_div);
 //   sorttable.makeSortable(table);
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