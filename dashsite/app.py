# -*- coding: utf-8 -*-
import base64

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
from components.algorithm_description import d
from pang.fileformat.json import reader as pangenomejson_reader
from pang.fileformat.json import writer as pangenomejson_writer


from networkx.readwrite import json_graph

app = dash.Dash(__name__)
app.title = 'pang'
app.layout = html.Div(
    children=[
        html.Div(
            id="header",
            children=[
                html.H1(
                    'Multiple sequence alignment analysis',
                    id='title',
                    className='ten columns',
                    style={'margin': '10px'}
                ),
                html.Img(
                    src=app.get_asset_url('logo.png'),
                    className='two columns',
                    style={'height':'60px',
                           'width':'60px',
                           'float':'right',
                           'margin': '10px'}
                )
            ],
            style={'backgroundColor': colors['light_shades'],
                   'height':'80px'},
            className='twelve columns row'
        ),
        html.Div(
            id="program",
            children=[
                        html.H3("Analysis", style={'margin-top':'0px'}),
                        html.Div(
                            children=[
                                    html.Div(
                                    id="program_run",
                                    children=[
                                        html.H5("Input files"),
                                        dcc.Upload(
                                            id="maf_upload",
                                            children=[
                                                html.Img(
                                                    src=app.get_asset_url('alignment.png'),
                                                    className='two columns file_upload_img',
                                                    style={'width': '50px', 'margin': '5px'}
                                                ),
                                                html.Div([html.A('Upload MAF file')],
                                                         className='ten columns'
                                                )
                                            ],
                                            className='file_upload'
                                        ),
                                        dcc.Upload(
                                            id="metadata_upload",
                                            children=[
                                                html.Img(
                                                    src=app.get_asset_url('information.png'),
                                                    className='two columns file_upload_img',
                                                    style={'width': '50px', 'margin': '5px'}
                                                ),
                                                html.Div([html.A('Upload metadata file')],
                                                         className='ten columns'
                                                )
                                            ],
                                            className='file_upload'
                                        ),
                                        html.H5("Algorithm version"),
                                        dcc.Dropdown(
                                            id='consensus_algorithm',
                                            options=[
                                                {'label': 'Simple', 'value': 'simple'},
                                                {'label': 'Tree', 'value': 'tree'},],
                                            value='tree',
                                            className='form_item'
                                        ),
                                        html.H5("Consensus generation options"),
                                        html.Div(
                                            id="consensus_params",
                                            children=[
                                                html.Div(
                                                    id='simple_consensus_params',
                                                    children=[
                                                        html.P("HBMIN value", className='form_item_label'),
                                                        dcc.Input(
                                                            id='hbmin',
                                                            placeholder='Enter hbmin value...',
                                                            min=0,
                                                            max=1.0,
                                                            type='number',
                                                            value=0.9
                                                        )],
                                                    style={'display': None},
                                                    className='form_item'
                                                ),
                                                html.Div(
                                                    id='tree_consensus_params',
                                                    children=[
                                                        html.P("Range", className='form_item_label'),
                                                        dcc.RangeSlider(
                                                            id='r',
                                                            marks={i: f'{i}' for i in range(0, 101, 5)},
                                                            min=0,
                                                            max=100,
                                                            value=[0, 100],
                                                            step=1,
                                                            className='form_item'
                                                        ),
                                                        html.Div(
                                                            children=[
                                                                html.P("Multiplier value", className='form_item_label'),
                                                                dcc.Input(
                                                                   id='multiplier',
                                                                   placeholder='Enter multiplier value...',
                                                                   min=0,
                                                                   type='number',
                                                                   value=0.7,
                                                                   className='form_item'
                                                                )],
                                                            style={'margin-top': '25px'}
                                                        ),
                                                        html.Div(
                                                            children=[
                                                                html.P("Stop value", className='form_item_label'),
                                                                dcc.Input(
                                                                    id='stop',
                                                                    placeholder='Enter stop value...',
                                                                    min=0,
                                                                    max=1.0,
                                                                    type='number',
                                                                    value=0.99,
                                                                className='form_item'
                                                                )]
                                                            ),
                                                            dcc.Checklist(
                                                                id='tree_consensus_options',
                                                                options=[
                                                                    {'label': 'Use re consensus', 'value': 're_consensus'},
                                                                    {'label': 'Anti fragmentatin node cutoff', 'value': 'anti_fragmentation'}
                                                                ],
                                                                values=['re_consensus', 'anti_fragmentation'],
                                                                className='form_item'
                                                            ),
                                                    ],
                                                    style={'display': None}
                                                )]
                                            ),
                                            html.H5("Additional options:"),
                                            dcc.Checklist(
                                                id='pang_options',
                                                options=[{'label': 'Generate fasta files', 'value': 'FASTA'}],
                                                values=[],
                                                                className='form_item'
                                            ),
                                            html.Div(
                                                [html.Button("Run pang",
                                                     id="pang_button",
                                                     className='button-primary form_item',
                                                ),
                                                html.A(html.Button("Download result as json",
                                                               id="json_download",
                                                               disabled=True,
                                                                   className='form_item'),
                                                    href='download_pangenome')],
                                                style={'margin': '10px'}
                                            ),
                                            html.Div(
                                                id='hidden_pang_result',
                                                style={'display': 'none'}
                                            )
                                        ],
                                    className='six columns'
                                        ),
                                        html.Div(
                                            id='program_info',
                                            children=[
                                                dcc.Markdown(d)
                                            ],
                                            className='six columns',
                                        )
                            ],
                        className='row'
                        )

            ],
            className='row'
        ),
        html.Div(
            id='visualisation',
            children=[
                html.H3("Visualisation"),
                html.Div(
                    id='data_upload',
                    children=[
                        dcc.Upload(
                                id='pangenome_json_upload',
                                children=html.Div([
                                    'Drag and Drop or ',
                                    html.A('Select Files')
                                ]),
                                style={
                                    # 'width': '100%',
                                    # 'height': '40px',
                                    'lineHeight': '40px',
                                    'borderWidth': '1px',
                                    'borderStyle': 'dashed',
                                    'borderRadius': '5px',
                                    'textAlign': 'center',
                                    'margin': '10px'
                                },
                                # Allow multiple files to be uploaded
                                multiple=True,
                                className='five columns'
                            ),
                        html.Div("or load example data: ",
                                 style={'textAlign':'center', 'lineHeight':'60px'},
                                 className='two columns'),
                        dcc.Dropdown(
                            options=[
                                {'label': 'Ebola', 'value': 'NYC'},
                                {'label': 'Mycoplasma', 'value': 'MTL'},
                                {'label': 'Proteins family', 'value': 'SF'}
                            ],
                            value='MTL',
                            style={
                                   'margin': '10px'},
                            className='five columns'
                        )
                    ],
                    style={}
                ),
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
            ],
            className='row'
        )
    ],
style={'margin':'0px',
       'padding':'0px'})


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
    anti_fragmentation_value = True if 'anti_fragmentation' in tree_consensus_options_values else False
    pangenome, json_path = run_pang(maf_contents,
                                    metadata_contents,
                                    fasta_option,
                                    consensus_algorithm_value,
                                    hbmin_value,
                                    r_value,
                                    multiplier_value,
                                    stop_value,
                                    re_consensus_value,
                                    anti_fragmentation_value)

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
