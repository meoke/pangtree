
consensuses_colors = ['#ff8000', '#602870', '#983352','yellow', 'purple', 'orange', 'brown']
document.addEventListener("DOMContentLoaded", function() {

   var cy = cytoscape({
    container: document.getElementById('cy'),
    elements: PoaGraph_elements,
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


    //poagraph_info
    var poagraph_info_div_id = document.getElementById("poagraph_info")
    poagraph_info_div_id.innerHTML += ('<p>' + 'Nodes count: ' + "\t\t" + poagraph.nodes_count + '<\p>' +
                                       '<p>' + 'Sources count:' + '\t\t' + poagraph.sources_count + '<\p>' +
                                       '<p>' + 'Consensuses count:' + '\t\t' + poagraph.consensuses_count + '<\p>' +
                                       '<p>' + 'Sequences per node (mean):' + '\t\t' + poagraph.mean_sources_per_node + '<\p>');

    //consensuses info
//  var divId = document.getElementById("sources_info")
  for(var i=0; i<consensuses.data.length; i++){
//    consensus_color = consensuses_colors[consensuses.data[i]['id']%consensuses_colors.length];
//    consensus_info_div = document.createElement('div')
//    consensus_info_div.setAttribute("class", "consensus_info")
//    consensus_info_div.innerHTML += ('<p style=\'color: ' + consensus_color +'\'>' + consensuses.data[i]['name'] + "</p>" +
//                        '<p>' + 'Nodes count:' + '\t' + consensuses.data[i]['nodes_count'] + '</p>');

//    var table = document.createElement('table');
//    table.className = "sortable";
//    var row = table.insertRow(0);
//    var head_id = row.insertCell(0);
//    var head_name = row.insertCell(1);
//    var head_title = row.insertCell(2);
//    var head_group = row.insertCell(3);
//    var head_bundle_id = row.insertCell(4);
//    var head_compatibility = row.insertCell(5);
//    head_id.innerHTML = "ID"
//    head_name.innerHTML = "Name"
//    head_title.innerHTML = "Title"
//    head_group.innerHTML = "Group"
//    head_bundle_id.innerHTML = "Bundle ID"
//    head_compatibility.innerHTML = "Compatibility"
//
//    for(var j=0;j<consensuses.data[i]['sources_map'].length;j++){
//        var row = table.insertRow(j+1);
//        var id = row.insertCell(0);
//        var name = row.insertCell(1);
//        var title = row.insertCell(2);
//        var group =row.insertCell(3)
//        var bundle_id = row.insertCell(4);
//        var compatibility = row.insertCell(5);
//        id.innerHTML = j
//        name.innerHTML = sources.data[j]['name']
//        title.innerHTML = sources.data[j]['title']
//        group.innerHTML = sources.data[j]['group_name']
//        bundle_id.innerHTML = sources.data[j]['bundle_ID']
//        compatibility.innerHTML = consensuses.data[i]['sources_map'][j]
//    }

//    compatibility_chart = get_density_graph(i);
//
//    consensus_info_div.appendChild(table)
//    consensus_info_div.appendChild(compatibility_chart)
//
//    divId.appendChild(consensus_info_div);
  }
});
