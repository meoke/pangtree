# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import shutil
import dash_table
import re

import pandas as pd
import flask
from app_style import colors, external_css
from pang_run import run_pang
from components import consensus_tree
from components import consensus_table
from pang.fileformat.json import reader as pangenomejson_reader
from pang.fileformat.json import writer as pangenomejson_writer



from networkx.readwrite import json_graph
df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/solar.csv')
app = dash.Dash(__name__)
app.title = 'pang'
app.layout = html.Div(
    children=[
        html.Div(
            children=[
                html.H2(
                    'Multiuliniowienie',
                    id='title'
                ),
                # html.Img(
                #     src="logo.png",
                # )
            ],
            style={'backgroundColor': colors['light_shades']},
        ),
        html.Div(
            id="program",
            children=[
                html.Div(
                    id="maf_upload_div",
                    children=[html.P("MAF file"),
                              dcc.Upload(
                                  id="maf_upload",
                                  children=[html.Div(
                                      children=[html.A("Select Files...")]
                                  )]
                              )]),
                html.Div(
                    id="metadata_upload_div",
                    children=[html.P("Metadata file"),
                              dcc.Upload(
                                  id="metadata_upload",
                                  children=[html.Div(
                                      children=[html.A("Select Files...")]
                                  )]
                              )]),
                dcc.Checklist(
                    id='pang_options',
                    options=[
                        {'label': 'Generate fasta files', 'value': 'FASTA'}
                    ],
                    values=[]
                ),
                dcc.Dropdown(
                    id='consensus_algorithm',
                    options=[
                        {'label': 'Simple', 'value': 'simple'},
                        {'label': 'Tree', 'value': 'tree'},
                    ],
                    value='tree'
                ),
                html.Div(
                    id="consensus_params",
                    children=[html.Div(
                        id='simple_consensus_params',
                        children=[html.P("HBMIN value"),
                                  dcc.Input(
                                      id='hbmin',
                                      placeholder='Enter hbmin value...',
                                      min=0,
                                      max=1.0,
                                      type='number',
                                      value=0.9
                                  )
                                  ],
                        style={'display': None}
                    ),
                        html.Div(
                            id='tree_consensus_params',
                            children=[html.P("Range"),
                                      dcc.RangeSlider(
                                          id='r',
                                          marks={i: f'{i}' for i in range(0, 101, 5)},
                                          min=0,
                                          max=100,
                                          value=[0, 100],
                                          step=1),
                                      html.Div(
                                          [html.P("Multiplier value"),
                                           dcc.Input(
                                               id='multiplier',
                                               placeholder='Enter multiplier value...',
                                               min=0,
                                               type='number',
                                               value=0.7
                                           )],
                            style={'margin-top': '25px'}
                                      ),
                                      html.Div(
                                          [html.P("Stop value"),
                                           dcc.Input(
                                               id='stop',
                                               placeholder='Enter stop value...',
                                               min=0,
                                               max=1.0,
                                               type='number',
                                               value=0.99
                                           )]
                                      ),
                                      dcc.Checklist(
                                          id='tree_consensus_options',
                                          options=[
                                              {'label': 'Use re consensus', 'value': 're_consensus'}
                                          ],
                                          values=['re_consensus']
                                      ),
                                      ],
                            style={'display': None}
                        )
                    ]
                ),
                html.Button("Run pang",
                            id="pang_button",
                            className='button-primary'
                            ),
                html.A(html.Button("Download result as json",
                                   id="json_download",
                                   disabled=True),
                       href='download_pangenome'),
                html.Div(
                    id='hidden_pang_result',
                    style={'display': 'none'}
                )
            ]
        ),
        html.Div(
            id='visualisation',
            children=[
                html.Div(
                    id='consensus_tree',
                    children=[
                        html.Div(
                            id='hidden_consensus_tree_data',
                            style={'display': 'none'}
                        ),
                        dcc.Graph(
                            id='consensus_tree_graph'
                        ),
                        dcc.Slider(
                            id='consensus_tree_slider',
                            min=0,
                            max=1,
                            marks={int(i) if i % 1 == 0 else i: '{}'.format(i) for i in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]},
                            step=0.01,
                            value=0.5,
                            dots=True
                        ),
                        html.Div(
                            id='consensus_tree_slider_value',
                            style={'margin-top': '25px'}
                        ),
                        html.Div(
                            id='hidden_consensuses_table_data',
                            style={'display': 'none'}
                        ),
                        html.Div(
                            dash_table.DataTable(
                                id='consensuses_table',
                                filtering=True,
                                sorting=True,
                                sorting_type="multi"
                            ),
                            style={'margin-top': '25px'},
                        )
                    ],
                    style={'display': 'none'}

                ),
                html.Div(
                    id='consensus_table',
                    children=[

                    ]
                ),
                html.Div(
                    id='blocks_graph'
                )
            ]
        )
    ])


@app.callback(
    dash.dependencies.Output('simple_consensus_params', 'style'),
    [dash.dependencies.Input('consensus_algorithm', 'value')])
def show_simple_consensus_algorithm_options(consensus_algorithm):
    if consensus_algorithm == 'simple':
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('tree_consensus_params', 'style'),
    [dash.dependencies.Input('consensus_algorithm', 'value')])
def show_tree_consensus_algorithm_options(consensus_algorithm):
    if consensus_algorithm == 'tree':
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('json_download', 'disabled'),
    [dash.dependencies.Input('hidden_pang_result', 'children')])
def download_json_button(_):
    return False


@app.callback(
    dash.dependencies.Output('hidden_pang_result', 'children'),
    [dash.dependencies.Input('pang_button', 'n_clicks')],
    [dash.dependencies.State('maf_upload', 'contents'),
     dash.dependencies.State('metadata_upload', 'contents'),
     dash.dependencies.State('pang_options', 'values'),
     dash.dependencies.State('consensus_algorithm', 'value'),
     dash.dependencies.State('hbmin', 'value'),
     dash.dependencies.State('r', 'value'),
     dash.dependencies.State('multiplier', 'value'),
     dash.dependencies.State('stop', 'value'),
     dash.dependencies.State('tree_consensus_options', 'values')
     ])
def call_pang(_,
              maf_contents,
              metadata_contents,
              pang_options_values,
              consensus_algorithm_value,
              hbmin_value,
              r_value,
              multiplier_value,
              stop_value,
              tree_consensus_options_values):
    fasta_option = True if 'FASTA' in pang_options_values else False
    re_consensus_value = True if 're_consensus' in tree_consensus_options_values else False
    pangenome, json_path = run_pang(maf_contents,
                                    metadata_contents,
                                    fasta_option,
                                    consensus_algorithm_value,
                                    hbmin_value,
                                    r_value,
                                    multiplier_value,
                                    stop_value,
                                    re_consensus_value)

    shutil.copy(json_path, "download/pangenome.json")
    return str(pangenomejson_writer.pangenome_to_json(pangenome))

@app.callback(
    dash.dependencies.Output('hidden_consensus_tree_data', 'children'),
    [dash.dependencies.Input('hidden_pang_result', 'children'),
     dash.dependencies.Input('consensus_tree_graph', 'clickData')],
    [dash.dependencies.State('hidden_consensus_tree_data', 'children')])
def update_consensus_graph_data(jsonified_pangenome, click_data, old_jsonfied_consensus_tree):
    old_tree = None
    if old_jsonfied_consensus_tree:
        old_tree = json_graph.tree_graph(eval(old_jsonfied_consensus_tree))

    print(click_data)
    jsonpangenome = pangenomejson_reader.json_to_jsonpangenome(jsonified_pangenome)
    tree = consensus_tree.get_tree(jsonpangenome, click_data, old_tree)
    jsonified_tree = json_graph.tree_data(tree, root=0)
    return str(jsonified_tree)


@app.callback(
    dash.dependencies.Output('consensus_tree', 'style'),
    [dash.dependencies.Input('consensus_tree_graph', 'figure')])
def show_graph(_):
    return {'display': 'block'}

@app.callback(
    dash.dependencies.Output('hidden_consensuses_table_data', 'children'),
    [dash.dependencies.Input('hidden_consensus_tree_data', 'children'),
     dash.dependencies.Input('consensus_tree_slider', 'value')],
    [dash.dependencies.State('hidden_pang_result', 'children')]
)
def update_consensuses_table(jsonified_consensus_tree, slider_value, jsonified_pangenome):
    tree = json_graph.tree_graph(eval(jsonified_consensus_tree))
    jsonpangenome = pangenomejson_reader.json_to_jsonpangenome(jsonified_pangenome)
    return str(consensus_table.get_consensus_table_data(jsonpangenome, tree, slider_value))


@app.callback(
    dash.dependencies.Output('consensuses_table', 'data'),
    [dash.dependencies.Input('hidden_consensuses_table_data', 'children')]
)
def update_consensuses_table_rows(jsonified_consensuses_table_data):
    a = re.sub("\s+", ",", jsonified_consensuses_table_data.strip())
    return consensus_table.get_table_rows(eval(a))

@app.callback(
    dash.dependencies.Output('consensuses_table', 'columns'),
    [dash.dependencies.Input('hidden_consensuses_table_data', 'children')]
)
def update_consensuses_table_columncs(jsonified_consensuses_table_data):
    a = re.sub("\s+", ",", jsonified_consensuses_table_data.strip())
    return consensus_table.get_table_columns(eval(a))

@app.callback(
    dash.dependencies.Output('consensus_tree_slider_value', 'children'),
    [dash.dependencies.Input('consensus_tree_slider', 'value')]
)
def show_slider_value(slider_value):
    return f"Current value: {slider_value}."


@app.callback(
    dash.dependencies.Output('consensus_tree_graph', 'figure'),
    [dash.dependencies.Input('hidden_consensus_tree_data', 'children'),
     dash.dependencies.Input('consensus_tree_slider', 'value')],
    [dash.dependencies.State('hidden_pang_result', 'children')])
def update_consensus_tree_graph(jsonified_consensus_tree, slider_value,  jsonified_pangenome):
    jsonpangenome = pangenomejson_reader.json_to_jsonpangenome(jsonified_pangenome)

    tree = json_graph.tree_graph(eval(jsonified_consensus_tree))
    return consensus_tree.get_consensus_tree_graph(jsonpangenome, tree, slider_value)


@app.server.route('/download_pangenome')
def download_csv():
    return flask.send_file('download/pangenome.json',
                           mimetype='text/csv',
                           attachment_filename='pangenome.json',
                           as_attachment=True)


for css in external_css:
    app.css.append_css({"external_url": css})

if __name__ == '__main__':
    app.run_server(debug=True)
