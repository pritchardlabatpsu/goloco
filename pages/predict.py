import dash
from dash import html, callback
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
#from dash_extensions.enrich import Dash, DashProxy, Output, Input, State, ServersideOutput, html, dcc, ServersideOutputTransform

dash.register_page(__name__)

SIDEBAR_STYLE = {
    "padding": "2rem 2rem",
    "background-color": "#f8f9fa",
}

INFERENCE_STYLE = {
    "background-color": "#F6F8FC",
    'borderWidth': '1px',
    'borderStyle': 'dashed',
    'borderRadius': '1px',
}

INFERENCE_STYLE_3 = {
    "color": "#f8f9fa",
}

INFERENCE_STYLE_2 = {
    "padding": "1rem 1rem",
    "background-color": "#F6F8FC",
    'borderWidth': '1px',
    'borderStyle': 'dashed',
    'borderRadius': '1px',
}

tab_style = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontSize': '20px',
    'fontWeight': 'bold',
    'font-family': 'Roboto, sans-serif',
    'width': '20%',
    'borderWidth': '1px',
    'borderStyle': 'solid',
    'borderRadius': '1px',
}

tab_selected_style = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': '#0f8de5',  #'#119DFF'
    'color': 'white',
    'font-family': 'Roboto, sans-serif',
    'fontSize': '20px',
    'fontWeight': 'bold',
    'padding': '6px',
    'width': '20%',
    'borderWidth': '1px',
    'borderStyle': 'solid',
    'borderRadius': '1px',
}

inference_frame = html.Div([

            html.P('''Use this tool to generate genome-wide loss-of-function predictions and data visualizations. Use option 1 to upload an file containing L200 CERES scores
            across your cell-line experiments, use option 2 to upload a csv of a prediction that was performed previously with option 1, or use option 3 to load a sample prediction
            if you are interested in testing the visualizations and predictions that this tool will generate.'''),

            html.Br(),

            dbc.Row([dbc.Col(

                dbc.CardGroup([

                dbc.Card([dbc.CardHeader(html.Div([dcc.Markdown('''###### Option 1: Submit a new job''')], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})),
                          dbc.CardBody([
                          html.Div([
                                       html.P('''Upload a csv or excel file of L200 CERES scores for your experiment to the box below. Select the run inference button to run your inference:                 '''),
                                        dcc.Upload(
                                        id='upload-L200-data',
                                        children=html.Div(['Drag and Drop or ',
                                                    html.A('Select L200 File')]),
                                                style={
                                                    'width': '100%',
                                                    'height': '40px',
                                                    'lineHeight': '40px',
                                                    'borderWidth': '1px',
                                                    'borderStyle': 'dashed',
                                                    'borderRadius': '5px',
                                                    'textAlign': 'center',
                                                    'margin': '0px',
                                                    'background-color': '#ffffff',
                                                },
                                                multiple=False
                                        ),

                                        dcc.Loading(id='file_upload_load', children = [html.Div([html.Label(id='file_upload_label', style={'color':'#548235'})])], type='default'),

                                        html.Div([
                                                    dbc.Button('Run Inference',
                                                    id='submit-inference',
                                                    color="primary",
                                                    disabled = True,
                                                    n_clicks = 0,
                                                    ),], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),
                                       ])])], style=INFERENCE_STYLE_2, outline=True),

                     dbc.Card([dbc.CardHeader(html.Div([dcc.Markdown('''###### Option 2: Upload a previous prediction''')], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})),
                               dbc.CardBody([
                               html.Div([
                                       html.P('''Upload a csv or excel file of a previous prediction to the box below. Select the load prediction button below to visualize your previous prediction data:    '''),
                                        dcc.Upload(
                                        id='upload-df_pred-data',
                                        children=html.Div(['Drag and Drop or ',
                                                html.A('Select Prediction File')]),
                                            style={
                                                'width': '100%',
                                                'height': '40px',
                                                'lineHeight': '40px',
                                                'borderWidth': '1px',
                                                'borderStyle': 'dashed',
                                                'borderRadius': '5px',
                                                'textAlign': 'center',
                                                'margin': '0px',
                                                'background-color': '#ffffff',
                                            },
                                            multiple=False
                                        ),
                                        dcc.Loading(id='file_upload_load', children = [html.Div(id='pred_upload_label', children=[html.Br()])], type='default'), #html.Div([html.Label(id='pred_upload_label', style={'color':'#548235'})])], type='default'

                                        html.Div([
                                                    dbc.Button('Load Prediction',
                                                    id='submit-prediction',
                                                    color="primary",
                                                    disabled = True,
                                                    n_clicks = None,
                                                    ),], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),

                                       ])])], style=INFERENCE_STYLE_2, outline=True),

                     dbc.Card([dbc.CardHeader(html.Div([dcc.Markdown('''###### Option 3: Run a sample visualization''')], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),),
                               dbc.CardBody([
                               html.Div([
                                       html.P('''Select the load sample button below to load visualizations of a sample prediction preloaded to our servers:                                                  '''),
                                       html.Br(),
                                       html.Div([
                                            dbc.Button('Load Sample',
                                            id='load-prediction',
                                            color="primary",
                                            n_clicks=0,),], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),
                                       html.Br(),
                                       html.Br(),
                                       html.Br(),
                                       ])])], style=INFERENCE_STYLE_2, outline=True),

                                    ], style=INFERENCE_STYLE),

                            ),
                     ]),

            html.Br(),

            dbc.Row([dbc.Col(), dbc.Col(html.Div([dbc.Progress(id='progress_bar', animated=True, striped=True, style={"visibility": "hidden", "height": "20px"}),], style={'height': '30px'}), width=6), dbc.Col(),]),

            dbc.Row([dbc.Col(), dbc.Col(html.Div([html.Label(id='eta_label', hidden=True),], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'},),), dbc.Col(),]),

            dbc.Row([dbc.Col(), 
                     dbc.Col(html.Div([
                        dbc.Button('Cancel',
                        id='cancel-inference',
                        color="primary",
                        n_clicks=0,),], style={'visibility': 'hidden', 'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),), 
                     dbc.Col()]),

                    ],style=SIDEBAR_STYLE),

gene_convert = html.Div([
                html.P('''Use this tool to convert your read-counts from CRISPR screens on L200 genes to CHRONOS scores.
                Converting raw reads to gene effect scores with the CHRONOS algorithm, developed by the 
                Broad Institute, requires at least 3 dataframes including a matrix of raw readcounts, a sequence 
                mapping, and an sgRNA guide mapping. Optionally, gene-level copy number calls can be submitted or 
                selected from the DepMap data to correct gene effect scores. Submit your data below in the correct 
                formats to convert to CHRONOS scores:'''),

                html.Br(),

            dbc.Row([dbc.Col(

                dbc.CardGroup([

                dbc.Card([dbc.CardHeader(html.Div([
                          dcc.Markdown('''###### 1) Readcounts Dataframe: ''')], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})),
                          dbc.CardBody([
                          html.Div([
                                       html.P("Upload a matrix as csv mapping count values to sgRNA sequences and sequence IDs with the following formats listed in the details below." , style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),
                                       dbc.Alert(
                                                [
                                                    html.Div([
                                                    html.I(className="bi bi-info-circle-fill me-2"),
                                                    html.Details([
                                                    html.Br(),
                                                    dbc.ListGroup(
                                                        [dbc.ListGroupItem([html.B('Indexes: '), 'unnamed index with records of sequenced entities including pDNA and biological replicates for your experiments']),
                                                         dbc.ListGroupItem([html.B('Columns: '), 'individual sequences of sgRNAs, first column in csv must be the unnamed index column']),
                                                         dbc.ListGroupItem([html.B('Values: '), 'number of reads counted for each sgRNA']),
                                                            ]),], open=False),
                                                    ],),
                                                ],
                                                color="info",
                                                className="d-flex align-items-center",
                                                style={'align-items': 'center'}
                                            ),
                                        dcc.Upload(
                                        id='upload-reads-data',
                                        children=html.Div(['Drag and Drop or ',
                                                    html.B('Select Readcounts File')]),
                                                style={
                                                    'width': '100%',
                                                    'height': '40px',
                                                    'lineHeight': '40px',
                                                    'borderWidth': '1px',
                                                    'borderStyle': 'dashed',
                                                    'borderRadius': '5px',
                                                    'textAlign': 'center',
                                                    'margin': '0px',
                                                    'background-color': '#ffffff',
                                                },
                                                multiple=False
                                        ),

                                        dcc.Loading(id='file_upload_load', children = [html.Div([html.Label(id='upload_reads_label', style={'color':'#548235'})], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})], type='default'),

                                       ])])], style=INFERENCE_STYLE_2, outline=True),

                     dbc.Card([dbc.CardHeader(html.Div([dcc.Markdown('''###### 2) Sequence Map Dataframe''')], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})),
                          dbc.CardBody([
                          html.Div([
                                       html.P("Upload a dataframe, as csv, mapping sequence IDs to cell lines, pDNA batches, and timepoints with the following four columns and precise headers as listed in the details below.", style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),
                                       dbc.Alert(
                                                [
                                                    html.Div([
                                                    html.I(className="bi bi-info-circle-fill me-2"),
                                                    html.Details([
                                                    html.Br(),
                                                    dbc.ListGroup(
                                                        [dbc.ListGroupItem([html.B('sequence_ID (str): '), 'sequenced entities, must match row indexes from readcounts matrix']),
                                                         dbc.ListGroupItem([html.B('cell_line_name (str): '), '"pDNA" for pDNA or cell-line name, each pDNA batch must have at least one pDNA measurement, biological replicates of the same cell-line should share the same cell_line_name']),
                                                         dbc.ListGroupItem([html.B('pDNA_batch (int or str): '), 'pDNA measurements in same batch are grouped and averaged, then used as reference for replicates assigned to batch']),
                                                         dbc.ListGroupItem([html.B('days (int): '), 'days post infection, ignored for pDNA']),
                                                            ],),
                                                        ], open= False, id='sequence-map-details')
                                                    ],),
                                                ],
                                                color="info",
                                                className="d-flex align-items-center",
                                                style={'align-items': 'center'}
                                            ),
                                        dcc.Upload(
                                        id='upload-sequence-map-data',
                                        children=html.Div(['Drag and Drop or ',
                                                    html.B('Select Sequencemap File')]),
                                                style={
                                                    'width': '100%',
                                                    'height': '40px',
                                                    'lineHeight': '40px',
                                                    'borderWidth': '1px',
                                                    'borderStyle': 'dashed',
                                                    'borderRadius': '5px',
                                                    'textAlign': 'center',
                                                    'margin': '0px',
                                                    'background-color': '#ffffff',
                                                    'display': 'flex',
                                                    'align-items': 'center',
                                                    'justify-content': 'center'
                                                },
                                                multiple=False
                                        ),

                                        dcc.Loading(id='file_upload_load', children = [html.Div([html.Label(id='upload_sequence_map_label', style={'color':'#548235'})], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})], type='default'),

                                       ])])], style=INFERENCE_STYLE_2, outline=True),

                     dbc.Card([dbc.CardHeader(html.Div([dcc.Markdown('''###### 3) Guide Map Dataframe''')], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})),
                          dbc.CardBody([
                          html.Div([
                                       html.P("Upload a dataframe mapping sgRNA sequences to genes as a csv with at least the two columns with precise headers as listed in the details below.", style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),
                                       dbc.Alert(
                                                [
                                                    html.Div([
                                                    html.I(className="bi bi-info-circle-fill me-2"),
                                                    html.Details([
                                                    html.Br(),
                                                    dbc.ListGroup(
                                                        [dbc.ListGroupItem([html.B('sgrna (str): '), 'sgRNA sequences, must match columns from readcounts matrix and can only appear once in this column']),
                                                         dbc.ListGroupItem([html.B('gene (str): '), 'the gene the sgRNA maps to in "geneSym (geneID)" format, i.e. "ACSL3 (2181)", should include the L200 genes']),
                                                            ],),
                                                        ], open= False, id='guide-map-details')
                                                    ],),
                                                ],
                                                color="info",
                                                className="d-flex align-items-center",
                                                style={'align-items': 'center'}
                                            ),
                                        dcc.Upload(
                                        id='upload-guide-map-data',
                                        children=html.Div(['Drag and Drop or ',
                                                    html.B('Select Guidemap File')]),
                                                style={
                                                    'width': '100%',
                                                    'height': '40px',
                                                    'lineHeight': '40px',
                                                    'borderWidth': '1px',
                                                    'borderStyle': 'dashed',
                                                    'borderRadius': '5px',
                                                    'textAlign': 'center',
                                                    'margin': '0px',
                                                    'background-color': '#ffffff',
                                                },
                                                multiple=False
                                        ),

                                        dcc.Loading(id='file_upload_load', children = [html.Div([html.Label(id='upload_guide_map_label', style={'color':'#548235'})])], type='default'),

                                       ])])], style=INFERENCE_STYLE_2, outline=True),
                                    ], style=INFERENCE_STYLE),

                            ),
                     ]),
            
            html.Br(),

            dcc.Loading(id='file_upload_load', children = [html.Div([html.Label(id='upload_all_convert_files', style={'color':'#548235'})], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})], type='default'),

            html.Div([
                      dbc.Button('Run CHRONOS Conversion',
                      id='submit-convert',
                      color="primary",
                      disabled=True,
                      n_clicks=0,),], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),
            
            html.Br(),    

            dbc.Row([dbc.Col(), dbc.Col(dcc.Loading(id = 'loading_convert', children = [html.Div(id='progress_bar_convert', style={'height': '30px'})], type='default'), width=6), dbc.Col(),]),

            dbc.Row([dbc.Col(), 
                     dbc.Col(html.Div([
                        dbc.Button('Cancel',
                        id='cancel-convert',
                        color="primary",
                        n_clicks=0,),], style={'visibility': 'hidden', 'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),), 
                     dbc.Col()]),
                ],style=SIDEBAR_STYLE),

layout = dbc.Container([
    html.Div([
    html.Br(),

    dbc.Row([dbc.Col(html.Div([
        dcc.Markdown('''#### Download Lossy Subsets'''),
        html.P("L100, L200, and L300 are compressed subsets of genes consisting of 100, 200, and 300 genes, respectively, which tunably predict thousands of other CRISPR gene effects. Genome-wide loss-of-function predictions can be made by measuring the gene effect of only these compressed subsets. You will find downloable links to the L100, L200, and L300 gene subsets below:"),
        html.A('L100_landmark_genes.csv', id='l100_download_link', href="#"),
        dcc.Download(id="l100_download_data"),
        html.Br(),
        html.A('L200_landmark_genes.csv', id='l200_download_link', href="#"),
        dcc.Download(id="l200_download_data"),
        html.Br(),
        html.A('L300_landmark_genes.csv', id='l300_download_link', href="#"),
        dcc.Download(id="l300_download_data"),
        html.Br(),],style=SIDEBAR_STYLE),),
             ]),
        
        html.Br(),

#    dbc.Row([dbc.Col(
        dcc.Tabs(id="tabs-styled-with-inline", value='infer-tab', children = [
            dcc.Tab(label='''Submit Inference''', value='infer-tab', style=tab_style, selected_style=tab_selected_style, children=inference_frame),
            dcc.Tab(label='''Gene Effect Calculator''', value='convert-tab', style=tab_style, selected_style=tab_selected_style, children=gene_convert),
        ]),

    #html.Div(id='tabs-content-inline', children = inference_frame),
    
    dbc.Modal(
        [
            dbc.ModalHeader("Warning"),
            dbc.ModalBody("Please upload your L200 CERES scores to the drag and drop box"),
            dbc.ModalFooter(
                dbc.Button(
                    "Close", id="close", color="primary", className="ms-auto", n_clicks=0
                )
            ),
        ],
        id="warning-popup",
        is_open=False,
    ),

    dbc.Modal(
        [
            dbc.ModalHeader("Warning"),
            dbc.ModalBody("Please upload your prediction csv prior to running your visualizations"),
            dbc.ModalFooter(
                dbc.Button(
                    "Close", id="close_2", color="primary", className="ms-auto", n_clicks=0
                )
            ),
        ],
        id="warning-popup-2",
        is_open=False,
    ),

    dbc.Modal(
        [
            dbc.ModalHeader("Information"),
            dbc.ModalBody("This inference may take over an hour. Please keep your browser and this page open as the inference runs. Once completed, you will be redirected to the overview page and a pop-up will prompt you to download your prediction. Select the button below to continue."),
            dbc.ModalFooter(
                dbc.Button(
                    "Continue", id="continue", color="primary", className="ms-auto", n_clicks=0
                )
            ),
        ],
        id="continue-popup",
        is_open=False,
    ),
    
    html.Br(),],     style={
        "transform": "scale(0.85)",
        "transform-origin": "top",
    },),

], fluid = True, style={'className': 'small'})

#---------------------------------------------Callbacks-----------------------------------------------------------------------------------------------------------------
@callback(
    Output("l100_download_data", "data"),
    Input("l100_download_link", "n_clicks"),
    prevent_initial_call=True
)
def return_l100(n_clicks):
    return dcc.send_file(
        "./depmap_files/L100_landmark_genes.csv"
        )

@callback(
    Output("l200_download_data", "data"),
    Input("l200_download_link", "n_clicks"),
    prevent_initial_call=True
)
def return_l200(n_clicks):
    return dcc.send_file(
        "./depmap_files/L200_landmark_genes.csv"
        )

@callback(
    Output("l300_download_data", "data"),
    Input("l300_download_link", "n_clicks"),
    prevent_initial_call=True
)
def return_l300(n_clicks):
    return dcc.send_file(
        "./depmap_files/L300_landmark_genes.csv"
        )

#@callback(Output('tabs-content-inline', 'children'),
#          Input('tabs-styled-with-inline', 'value'),
#          prevent_initial_call=True)
#def render_content(tab):
##    if tab == 'infer-tab':
 #       return inference_frame
 #   elif tab == 'convert-tab':
 #       return gene_convert