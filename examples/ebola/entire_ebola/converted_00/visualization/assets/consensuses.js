
consensuses_colors = ['#ff8000', '#602870', '#983352','yellow', 'purple', 'orange', 'brown']
document.addEventListener("DOMContentLoaded", function() {
    treshold_index=0;
    //poagraph_info
    var poagraph_info_div_id = document.getElementById("poagraph_info")
    poagraph_info_div_id.innerHTML += ('<p>' + 'Nodes count: ' + "\t\t" + poagraph.nodes_count + '<\p>' +
                                       '<p>' + 'Sources count:' + '\t\t' + poagraph.sources_count + '<\p>' +
                                       '<p>' + 'Sequences per node (mean):' + '\t\t' + poagraph.mean_sources_per_node + '<\p>');

    //consensuses info
  var divId = document.getElementById("sources_info")
  //for(var i=0; i<consensuses.data.length; i++){
    //consensus_color = consensuses_colors[consensuses.data[i]['id']%consensuses_colors.length];
    consensus_info_div = document.createElement('div')
    consensus_info_div.setAttribute("class", "consensus_info")
    //consensus_info_div.innerHTML += ('<p style=\'color: ' + consensus_color +'\'>' + consensuses.data[i]['name'] + "</p>" +
    //                    '<p>' + 'Nodes count:' + '\t' + consensuses.data[i]['nodes_count'] + '</p>');

    var table = document.createElement('table');
    table.className = "sortable";
    var row = table.insertRow(0);
    var head_id = row.insertCell(0);
    var head_name = row.insertCell(1);
    var head_title = row.insertCell(2);
    var head_group = row.insertCell(3);
    var head_bundle_id = row.insertCell(4);
    var consensuses_count = consensuses.c[treshold_index].data.length;
    head_id.innerHTML = "ID"
    head_name.innerHTML = "Name"
    head_title.innerHTML = "Title"
    head_group.innerHTML = "Group"
    head_bundle_id.innerHTML = "Bundle ID"

    for(var i=0; i< consensuses_count; i++){
    var cons_compatibility = row.insertCell(5+i);
    cons_compatibility.innerHTML = consensuses.c[treshold_index].data[i].name;
    }

    for(var j=0;j<sources.data.length;j++){
        var row = table.insertRow(j+1);
        var id = row.insertCell(0);
        var name = row.insertCell(1);
        var title = row.insertCell(2);
        var group =row.insertCell(3)
        var bundle_id = row.insertCell(4);
        id.innerHTML = j
        name.innerHTML = sources.data[j]['name']
        title.innerHTML = sources.data[j]['title']
        group.innerHTML = sources.data[j]['group_name']

        b_ID = -1
        cons_max = -1
        for(var k=0;k<consensuses.c[treshold_index].data.length; k++){
            var compatibility = row.insertCell(5+k);
            c = consensuses.c[treshold_index].data[k]['sources_compatibility'][j]
            compatibility.innerHTML = c;
            if(c > cons_max){
                b_ID = k
                cons_max = c
            }
        }
        bundle_id.innerHTML = b_ID

    }
//    for(var j=0;j<consensuses.data[i]['sources_compatibility'].length;j++){
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
//        compatibility.innerHTML = consensuses.data[i]['sources_compatibility'][j]
//    }

    //compatibility_chart = get_density_graph(i);

    consensus_info_div.appendChild(table);
    //consensus_info_div.appendChild(compatibility_chart)

    divId.appendChild(consensus_info_div);
  //}
});
