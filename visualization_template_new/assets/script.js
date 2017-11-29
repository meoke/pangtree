consensuses_colors = ['#ff8000', '#602870', '#983352','yellow', 'purple', 'orange', 'brown']

document.addEventListener("DOMContentLoaded", function() {
    add_poagraph_info();
    add_consensuses_options();
//    draw_visualization();
//    show_consensuses_tree();
//    show_sources_info_table();


//    //poagraph_info
//    var poagraph_info_div_id = document.getElementById("poagraph_info")
//    poagraph_info_div_id.innerHTML += ('<p>' + 'Nodes count: ' + "\t\t" + poagraph.nodes_count + '<\p>' +
//                                       '<p>' + 'Sources count:' + '\t\t' + poagraph.sources_count + '<\p>' +
//                                       '<p>' + 'Consensuses count:' + '\t\t' + poagraph.consensuses_count + '<\p>' +
//                                       '<p>' + 'Sequences per node (mean):' + '\t\t' + poagraph.mean_sources_per_node + '<\p>');
//
//    //consensuses info
//  var divId = document.getElementById("sources_info")
//  for(var i=0; i<consensuses.data.length; i++){
//    consensus_color = consensuses_colors[consensuses.data[i]['id']%consensuses_colors.length];
//    consensus_info_div = document.createElement('div')
//    consensus_info_div.setAttribute("class", "consensus_info")
//    consensus_info_div.innerHTML += ('<p style=\'color: ' + consensus_color +'\'>' + consensuses.data[i]['name'] + "</p>" +
//                        '<p>' + 'Nodes count:' + '\t' + consensuses.data[i]['nodes_count'] + '</p>');
//
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
//
//    compatibility_chart = get_density_graph(i);
//
//    consensus_info_div.appendChild(table)
//    consensus_info_div.appendChild(compatibility_chart)
//
//    divId.appendChild(consensus_info_div);
//  }

});



function add_poagraph_info(){
    name_p = document.getElementById('name');
    name_p.innerHTML += info.name;

    consensus_algorithm_p = document.getElementById('consensus_algorithm');
    consensus_algorithm_p.innerHTML += info.consensus_algorithm;

    running_time_p = document.getElementById('running_time');
    running_time_p.innerHTML += info.running_time;

    nodes_count_p = document.getElementById('nodes_count');
    nodes_count_p.innerHTML += info.nodes_count;

    sequences_per_node_p = document.getElementById('sequences_per_node');
    sequences_per_node_p.innerHTML += info.sequences_per_node;
}

function add_consensuses_options(){
    for(i=0; i<info.levels.length;i++)
    {
        let level_value = info.levels[i];
        var radioBtn_p = $('<p></p>')

        var radioBtn = $('<input type="radio" name="level" value="' + String(level_value) + '">');
        radioBtn.click( function(){
            show_sources_info_table(level_value);
        });
        radioBtn.appendTo(radioBtn_p);

        var radioBtn_label = $('<label>' + String(level_value) + '</label>')
        radioBtn_label.appendTo(radioBtn_p);

        radioBtn_p.appendTo('#consensuses_info');
    }

}

function draw_visualization(){
    //alert('vis');
}

function show_consensuses_tree(){
    //alert('tree');
}

function show_sources_info_table(level_value){
    var divId = document.getElementById("sources_info")
}