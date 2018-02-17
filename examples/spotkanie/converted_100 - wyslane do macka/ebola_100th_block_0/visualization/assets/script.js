consensuses_colors = ['yellow', 'green', 'orange', 'red', 'blue', 'violet', 'black', 'grey', 'pink']

document.addEventListener("DOMContentLoaded", function() {
    add_poagraph_info();
    add_consensuses_options();
    draw_visualization();
    //show_consensuses_tree();
    //show_sources_info_table();

});



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
    for(i=0; i<info.levels.length;i++)
    {
        let level_value = info.levels[i];
        let level_index = i;
        var radioBtn_p = $('<p></p>')

        var radioBtn = $('<input type="radio" name="level" value="' + String(level_value) + '">');
        radioBtn.click( function(){
            show_sources_info_table(level_index);
        });
        radioBtn.appendTo(radioBtn_p);

        var radioBtn_label = $('<label>' + String(level_value) + '</label>')
        radioBtn_label.appendTo(radioBtn_p);

        radioBtn_p.appendTo('#consensuses_info');
    }

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
          'curve-style': 'bezier',
          'label': 'data(title)',
          'line-color': get_edge_label_color, //to differentiate consensuses
          'width': get_edge_width, //to show popular edges
      }
    },{
      selector: '.aligned',
      style: {
        'line-style': 'dotted',
        'curve-style': 'haystack'
      }
    },{
      selector: '.consensus',
      style: {
        'curve-style': 'bezier'
      }
    }
    ],
    layout: {
      name: 'preset'
    }
  });

  function get_edge_label_color(ele) {

    consensus_id = ele.data('consensus');
    //console.log(consensus_id);
    if(consensus_id === -1){
      return 'white';
    }
    c = consensuses_colors[consensus_id%consensuses_colors.length];
    return c;
  }

  function get_edge_width(ele) {
    edge_width = Math.abs(Math.log(ele.data('weight') * 50))
    return edge_width.toString() + 'px';
  }

  function get_node_width(ele) {
    node_width = Math.abs(Math.log(ele.data('weight') * 1000))
    return node_width.toString() + 'px';
  }

}

function show_consensuses_tree(){
//    var example_tree = "(((EELA:0.150276,CONGERA:0.213019):0.230956,(EELB:0.263487,CONGERB:0.202633):0.246917):0.094785,((CAVEFISH:0.451027,(GOLDFISH:0.340495,ZEBRAFISH:0.390163):0.220565):0.067778,((((((NSAM:0.008113,NARG:0.014065):0.052991,SPUN:0.061003,(SMIC:0.027806,SDIA:0.015298,SXAN:0.046873):0.046977):0.009822,(NAUR:0.081298,(SSPI:0.023876,STIE:0.013652):0.058179):0.091775):0.073346,(MVIO:0.012271,MBER:0.039798):0.178835):0.147992,((BFNKILLIFISH:0.317455,(ONIL:0.029217,XCAU:0.084388):0.201166):0.055908,THORNYHEAD:0.252481):0.061905):0.157214,LAMPFISH:0.717196,((SCABBARDA:0.189684,SCABBARDB:0.362015):0.282263,((VIPERFISH:0.318217,BLACKDRAGON:0.109912):0.123642,LOOSEJAW:0.397100):0.287152):0.140663):0.206729):0.222485,(COELACANTH:0.558103,((CLAWEDFROG:0.441842,SALAMANDER:0.299607):0.135307,((CHAMELEON:0.771665,((PIGEON:0.150909,CHICKEN:0.172733):0.082163,ZEBRAFINCH:0.099172):0.272338):0.014055,((BOVINE:0.167569,DOLPHIN:0.157450):0.104783,ELEPHANT:0.166557):0.367205):0.050892):0.114731):0.295021)"
    var tree = d3.layout.phylotree().svg(d3.select("#tree_display"));
    tree(d3.layout.newick_parser(newick)).layout();
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
    var consensuses_count = consensuses.c[level_index].data.length;
    head_id.innerHTML = "ID"
    head_name.innerHTML = "Name"
    head_title.innerHTML = "Title"
    head_group.innerHTML = "Group"
    head_bundle_id.innerHTML = "Bundle ID"

    for(var i=0; i< consensuses_count; i++){
    var cons_compatibility = header_row.insertCell(5+i);
    cons_compatibility.innerHTML = consensuses.c[level_index].data[i].name;
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
        bundle_id.innerHTML = sources.data[j]['bundle_ID'][level_index];


        for(var k=0;k<consensuses.c[level_index].data.length; k++){
            var compatibility = row.insertCell(5+k);
            c = consensuses.c[level_index].data[k]['sources_compatibility'][j]
            compatibility.innerHTML = c;
        }
    }
 table.appendChild(body)
 table.createTFoot()
    consensus_info_div.appendChild(table)
    sources_info_div.appendChild(consensus_info_div);
    sorttable.makeSortable(table);
}