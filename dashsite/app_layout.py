import dash_html_components as html
import dash_core_components as dcc
import dash_table

from app_style import colors
from components.algorithm_description import general_info, table_info_markdown


def get_layout(get_url_function):
    return html.Div(
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
                    src=get_url_function('logo.png'),
                    className='two columns',
                    style={'height': '60px',
                           'width': '60px',
                           'float': 'right',
                           'margin': '10px'}
                )
            ],
            style={'backgroundColor': colors['accent'],
                   'height':'80px',
                   'width': '100%',
                   'color': colors['text']},
            className='twelve columns row'
        ),
        html.Div(
            id='body',
            children=[
                html.Div(
                    id="program",
                    children=[
                        # html.H3("Analysis", style={'margin-top': '0px'}),
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
                                                    src=get_url_function('alignment.png'),
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
                                                    src=get_url_function('information.png'),
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
                                            id='algorithm',
                                            options=[
                                                {'label': 'Simple', 'value': 'simple'},
                                                {'label': 'Tree', 'value': 'tree'}, ],
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
                                                                    value=1,
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
                                                                {'label': 'Anti fragmentatin node cutoff',
                                                                 'value': 'anti_fragmentation'}
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
                                        ),
                                        html.Div(
                                            id='hidden_last_clicked',
                                            style={'display': 'none'}
                                        )
                                    ],
                                    className='six columns'
                                ),
                                html.Div(
                                    id='program_info',
                                    children=[
                                        dcc.Markdown(general_info)
                                    ],
                                    className='six columns',
                                    style={'backgroundColor': colors['cold_background']}
                                )
                            ],
                            className='row',
                            style={'padding': '100px 0px 0px 0px'}
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
                                    id='pangenome_upload',
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
                                         style={'textAlign': 'center', 'lineHeight': '60px'},
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
                                ),
                                html.Button("Load pangenome",
                                            id="load_pangenome",
                                            className='button-primary form_item',
                                            )
                            ],
                            style={}
                        ),
                        html.Div(
                            id='program_parameters',
                            children=[
                                html.Div(
                                    id='hidden_program_parameters',
                                    style={'display': 'none'}
                                ),
                                html.Div(
                                    id='program_parameters_display'
                                )
                            ],
                            style={'width': '100%'}
                        ),
                        html.Div(
                            id='consensus_tree',
                            children=[
                                html.Div(
                                    id='hidden_consensus_tree_data',
                                    style={'display': 'none'}
                                ),
                                html.Div(
                                    id='tree',
                                    children=[
                                        html.Div(
                                            id='graphics',
                                            children=[
                                                dcc.Graph(
                                                    id='consensus_tree_graph',
                                                    style={'height': '1000px', 'width': 'auto'}
                                                ),
                                                html.Div(
                                                    [html.Div(
                                                        dcc.Slider(
                                                            id='consensus_tree_slider',
                                                            min=0,
                                                            max=1,
                                                            marks={int(i) if i % 1 == 0 else i: '{}'.format(i) for i in
                                                                   [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]},
                                                            step=0.01,
                                                            value=0.5,
                                                            dots=True
                                                        ),
                                                        style={'margin-top': '1%'},
                                                        className='ten columns'
                                                    ),
                                                        html.P(
                                                            id='consensus_tree_slider_value',
                                                            style={'font-size': 'large'},
                                                            className='two columns'
                                                        )],
                                                    className='row',
                                                    style={'margin-left': '3%',
                                                           'margin-right': '2%',
                                                           'margin-top': '-7%'}
                                                ),
                                            ],
                                            className='nine columns'
                                        ),
                                        html.Div(
                                            id='tree_info',
                                            children=[
                                                html.Div(
                                                    "<Compatibility charts placeholder>"
                                                ),
                                                html.H5(
                                                    id='consensus_node_details_header'
                                                ),
                                                dash_table.DataTable(
                                                    id='consensus_node_details',
                                                    columns=[{"name": i, "id": i} for i in ["ID", "Name", "Group"]],
                                                    style_table={
                                                        'maxHeight': '800',
                                                        'overflowY': 'scroll'
                                                    },
                                                    style_cell={'textAlign': 'left'},
                                                    sorting=True
                                                )
                                            ],
                                            style={'padding-top': '7%', 'padding-right': '2%'},
                                            className='three columns'
                                        )
                                    ],
                                    className='row'
                                ),
                                html.Div(
                                    id='hidden_consensuses_table_data',
                                    style={'display': 'none'}
                                ),
                                html.Div(
                                    children=[html.Div(
                                        id='consensus_table_info',
                                        children=[
                                            dcc.Markdown(table_info_markdown, className='ten columns'),
                                            html.A(html.Button("Download table as csv",
                                                               id="csv_download",
                                                               disabled=False,
                                                               className='form_item two columns'),
                                                   href='download_csv'),
                                            html.Div(id='hidden_csv_generated',
                                                     style={'display': 'none'})
                                        ],
                                        style={'padding': '2%'}
                                    ),
                                        dash_table.DataTable(
                                            id='consensuses_table',
                                            sorting=True,
                                            sorting_type="multi"
                                        )],
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
                        ),
                        html.Div(
                            id='pangraph',
                            children=[
                                html.Div(
                                    id='hidden_pangraph_points',
                                    style={'display': 'none'}
                                ),
                                html.Div(
                                    id='hidden_pangraph_traces',
                                    style={'display': 'none'}
                                ),
                                html.Div(
                                    id='pangraph_display',
                                    children=[
                                        dcc.Graph(
                                            id='pangraph_graph',
                                            style={'height': '300px', 'width': 'auto'}
                                        )
                                    ],
                                    style={'display': 'none'}
                                )
                            ]
                        )
                    ],
                    className='row'
                )
            ],
            style={'padding': '0px 30px 0px 30px'}
        )

    ],
    style={'padding': '0px',
            'margin': '0px'})
