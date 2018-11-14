# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
colors = {
    'light_shades': '#F2F0F2',
    'light_accent': '#79A6A4',
    'main_brand': '#73C700',
    'dark_accent': '#7D603E',
    'dark_dashes': '#2B2730'
}

app.layout = html.Div(style={'backgroundColor': colors['light_shades']}, children=[
    html.Div([
            html.H2(
                'Multiuliniowienie',
                id='title'
            ),
            html.Img(
                src="logo.png"
            )
        ],
            # className="banner",
            style={
            # 'textAlign': 'center',
            'backgroundColor': colors['main_brand']
        }
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
                id='opt',
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
                id="consensus_options"
            ),
            html.Button("Run pang",
                        id="pang_button"),
            html.Div(
                id='program_result'
            )
            ]),
    html.Div(
        id='visualisation'
    )
])


@app.callback(
    dash.dependencies.Output('consensus_options', 'children'),
    [dash.dependencies.Input('consensus_algorithm', 'value')])
def show_consensus_algorithm_options(consensus_algorithm):
    if consensus_algorithm == 'simple':
        return [html.P("HBMIN value"),
                dcc.Input(
                placeholder='Enter hbmin value...',
                min=0,
                max=1.0,
                type='number',
                value=0.9
                )
                ]
    elif consensus_algorithm == 'tree':
        return [html.P("Range"),
                dcc.RangeSlider(
                    marks={i: f'{i}' for i in range(0, 101, 5)},
                    min=0,
                    max=100,
                    value=[0, 100],
                    step=1),
                html.Div(
                    [html.P("Multiplier value"),
                     dcc.Input(
                     placeholder='Enter multiplier value...',
                     min=0,
                     type='number',
                     value=0.7
                     )]
                ),
                html.Div(
                    [html.P("Stop value"),
                     dcc.Input(
                        placeholder='Enter stop value...',
                        min=0,
                        max=1.0,
                        type='number',
                        value=0.99
                    )]
                ),
                dcc.Checklist(
                    options=[
                        {'label': 'Use re consensus', 'value': 're_consensus'}
                    ],
                    values=['re_consensus']
                ),
                ]

# @app.callback(
#     dash.dependencies.Output('program_result', 'children'),
#     [(dash.dependencies.Input('pang_button', 'n_clicks'))],
#     [dash.dependencies.State('opt', 'values')])
# def call_pan(n_clicks, values):
#     #maf file
#     #metadata file
#     #fasta
#     #vis
#     #consensus
#     #hbmin
#     #r
#     #multiplier
#     #stop
#     #re consensus
#     return 'The program value was "{}" and the button has been clicked {} times'.format(
#         values,
#         n_clicks
#     )

def parse_contents(content, filename):
    return html.Div(f"File: {content}, {filename}")

@app.callback(
    dash.dependencies.Output('program_result', 'children'),
    [dash.dependencies.Input('pang_button', 'n_clicks')],
     [dash.dependencies.State('maf_upload', 'contents'),
     dash.dependencies.State('metadata_upload', 'contents'),
     dash.dependencies.State('opt', 'values'),
     dash.dependencies.State('consensus_algorithm', 'value')])
def update_output(n_clicks, value, ff):
    return html.Div(f"F{n_clicks}, c{value}, f{ff}")
#
# @app.callback(dash.dependencies.State('test', 'children'),
#               [dash.dependencies.State('maf_upload', 'contents')],
#               [dash.dependencies.Input('maf_upload', 'contents')],
#               [dash.dependencies.State('maf_upload', 'filename')])
# def update_output(list_of_contents, list_of_names):
#     print("AA", list_of_contents)
#     print("BB", list_of_names)
#     return html.Div(f"{list_of_names}")
#     # if list_of_contents is not None:
#     #     children = [
#     #         parse_contents(c, n) for c, n in
#     #         zip(list_of_contents, list_of_names)]
#     #     return children


external_css = [
    "https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css"
]

for css in external_css:
    app.css.append_css({"external_url": css})

if __name__ == '__main__':
    app.run_server(debug=True)
