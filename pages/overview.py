import dash
from dash import html, callback
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash import dash_table
from dash.dash_table.Format import Format
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import matplotlib

full_color_list = mcolors.CSS4_COLORS
colormap_list = matplotlib.pyplot.colormaps()

dash.register_page(__name__)

SIDEBAR_STYLE = {
    "padding": "1rem 1rem",
    "background-color": "#f8f9fa",
}

tab_style = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontSize': '20px',
    'font-family': "Helvetica",
    'width': '20%',
    'borderWidth': '1px',
    'borderStyle': 'solid',
    'borderRadius': '1px',
}

tab_selected_style = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': '#0f8de5',
    'color': 'white',
    'font-family': 'Helvetica',
    'fontSize': '20px',
    'padding': '6px',
    'width': '20%',
    'borderWidth': '1px',
    'borderStyle': 'solid',
    'borderRadius': '1px',
}


df_dummy_1 = pd.DataFrame(columns=['Gene', 'Gene Category', 'Average', 'Standard Deviation', 'CERES Pred', 'Z-Score Pred'])

table_frame = html.Div([
        html.P('Select experiments below to generate a table of predictions',
        style={'font-size': '1.2rem', 'font-style': 'italic'}),
        html.Div([dcc.Dropdown(id='choose_experiment_1', multi=True, placeholder='Select experiments...')], style={"width": "40%"}),
        html.Br(),
        dcc.Loading(id = 'table_loading', children =[
        html.Div(id='output-inference-table', children = dash_table.DataTable(data=df_dummy_1.to_dict('records'),
                                              columns=[{'name': 'Gene', 'id': 'Gene', 'type':'text'},
                                                       {'name': 'Gene Category', 'id': 'Gene Category', 'type': 'text'},
                                                       {'name': 'Average', 'id': 'Average', 'type': 'numeric', 'format':Format(precision=2)},
                                                       {'name': 'Standard Deviation', 'id': 'Standard Deviation', 'type': 'numeric', 'format':Format(precision=2)},
                                                       {'name': 'CERES Prediction', 'id': 'CERES Pred', 'type': 'numeric', 'format':Format(precision=2)},
                                                       {'name': 'Z Score', 'id': 'Z-Score Pred', 'type': 'numeric', 'format':Format(precision=2)}],
                                              filter_action='native',
                                              sort_action='native',
                                              sort_mode='multi',
                                              style_table={
                                                    'height': 800,
                                                    'overflow': 'hidden',
                                                    'textOverflow': 'ellipsis',
                                              },
                                              style_header={
                                                            'backgroundColor': 'dodgerblue',
                                                            'font_family': 'sans-serif',
                                                            'font_size': '18px',
                                                            'paddingLeft': '20px',
                                                            'padding': '5px',
                                                            'color': 'white',
                                                            'fontWeight': 'bold'
                                                        },
                                              style_cell={'textAlign': 'left',                               
                                                          'font_family': 'Arial',
                                                          'font_size': '16px',
                                                          'paddingLeft': '20px',
                                                          'padding': '5px'},
                                              style_data={
                                                    'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                                                    'overflow': 'hidden',
                                                    'tex4tOverflow': 'ellipsis',
                                                    'minWidth': 50,
                                              }
                                              ))], type = 'default'),
        ],style=SIDEBAR_STYLE)

distributions_frame = html.Div([
    dbc.Row([
        dbc.Col(
            html.Div([
            html.P('Select an experiment below to visualize the distribution of CERES scores and z-scores of predictions',
            style={'font-size': '1.2rem', 'font-style': 'italic'}),
            html.Div([dcc.Dropdown(id='choose_experiment_2', placeholder='Select an experiment...'),
                      html.Br(),
                      dcc.Dropdown(id='distribution_color', placeholder='Select a colormap...', options = [{'label':k, 'value':k} for k in colormap_list])], style={"width": "40%"}),
            html.Br(),
            dbc.Checklist(
                        options = [{"label":'Conditional Essential', "value": "conditional essential"},
                        {"label":'Nonessential', "value": "common nonessential"},
                        {"label":'Essential', "value": "common essential"},
                        ], 
                        inline=True, 
                        inputStyle={"margin-right": "6px", "margin-left": "15px"},
                        value=[],
                        id='gene_cat_checklist_1'),

            html.Br(),

            dbc.Card([dbc.FormGroup([
            dbc.Row([
                dbc.Col(
                        html.Div([
                            dcc.Loading(id='ceres_distribution', children =[dcc.Graph(id='total_dist_graph')], type='default'),],
                                )
                    ),
                dbc.Col(html.Div([
                        dcc.Loading(id='ceres_distribution', children=[dcc.Graph(id='z_score_dist_graph')], type='default'),],
                                )
                        ),

                ]),

            html.Br(),

            dbc.Row([
                dbc.Col(
                        html.Div([
                            dcc.Loading(id='pdf_distribution', children=[dcc.Graph(id='total_dist_graph_pdf')], type='default'),],
                                )
                    ),
                dbc.Col(html.Div([
                        dcc.Loading(id='pdf_distribution_zscore', children=[dcc.Graph(id='z_score_dist_graph_pdf')], type='default'),],
                                )
                        ),

                ]),
            ])], style={'padding': '10px'}),

        ],),
    ),
        dbc.Col(
            dbc.Card([dbc.FormGroup([
            html.Div([
                    html.H4("Explore Genes", className="display-5"),
                    html.P('Click on a bin from the histograms on the left to view genes within that bin:', style={'font-style': 'italic'}),
                    html.Div(id='hist-gene-table',
                             children = [dcc.Loading(id='load_table', children=[dbc.Table.from_dataframe(df_dummy_1[['Gene', 'CERES Pred', 'Z-Score Pred']], striped=True, bordered=True, hover=True, size='sm')], type='default'),], 
                                 style={
                                            "height": "1000px",  # set the height of the div
                                            "overflow": "scroll",  # enable scrolling
                                        },),
                            ],
                         )])], style={'padding': '10px'}), width=3),

    ]),], style=SIDEBAR_STYLE)

z_score_frame = dbc.Row([dbc.Col(html.Div([
        html.P('Select an experiment below to visualize predictions of selective essential genes with gene effect score and z-score thresholds', 
        style={'font-size': '1.2rem', 'font-style': 'italic'}),

        dbc.Row([
            dbc.Col(html.Div([
                    dbc.Card([dbc.FormGroup([html.H6('Primary Histograms:', className='card-title'),
                        dbc.Label('Select Experiment:'),
                        dbc.Select(id='choose_experiment_3'),
                        html.Br(),
                        dbc.Label('Select Experiment Color:'),
                        dbc.Select(id='experiment_c',
                                   options = [{'label':k, 'value':k} for k in full_color_list],
                                   value='red'),
                        html.Br(),
                        dbc.Label('Select Primary Histogram Color:'),
                        dbc.Select(id='primary_c',
                                   options = [{'label':k, 'value':k} for k in full_color_list],
                                   value='green'),
                        html.Br(),
                        html.Br(),
                        dbc.Label('Select Max Genes:'),
                        dbc.Select(
                        id="select_number_genes",
                            options=[{'label':'  10', 'value': 10},
                                    {'label':'  20', 'value': 20},
                                    {'label':'  50', 'value': 50},
                                    {'label':'  100', 'value': 100},
                                    {'label':'  200', 'value': 200},
                                    {'label':'  500', 'value': 500},
                                    {'label':'  1000', 'value': 1000},
                                    {'label':'  2000', 'value': 2000},],
                            value=100,
                        ),
                        #html.Br(),
                        html.Br(),
                        dbc.Label('Select CERES Threshold:'),
                        dcc.Slider(-2, 1, 0.1, value=-0.7, id='sig_genes_threshold_1',
                        marks=dict(zip([i for i in np.arange(-2.5, 1.5, 0.5)], ["{:.1f}".format(i) for i in np.arange(-2.5, 1.5, 0.5)])),
                        tooltip={"placement": "bottom", "always_visible": False}),
                        #html.Br(),
                        dbc.Label('Select Average Range Across Cell-lines:'),
                        dcc.RangeSlider(-2, 1, 0.1, value=[-1, -0.3], id='z-score_threshold_1',
                        marks=dict(zip([i for i in np.arange(-2.5, 1.5, 0.5)], ["{:.1f}".format(i) for i in np.arange(-2.5, 1.5, 0.5)])),
                        tooltip={"placement": "bottom", "always_visible": False}
                        ),
                        html.Br(),
                        html.H6('Secondary Histograms:', className= 'card-title'),
                        dbc.Label('Select Category:'),
                        dbc.Select(
                                    id="select_primary_cat",
                                    options=[{'label':'  Lineage', 'value':'lineage'},
                                            {'label':'  Lineage Subtype', 'value':'lineage_subtype'},
                                            {'label':'  Collection Site', 'value':'sample_collection_site'},
                                            {'label':'  Disease', 'value':'disease'},
                                            {'label':'  Disease Subtype', 'value':'disease_subtype'},
                                            {'label':'  Origin', 'value':'primary_or_metastasis'},],
                                    ),
                        html.Br(),
                        dbc.Label('Select Subcategory:'),
                        dbc.Select(id="select_sub_cat",),
                        html.Br(),
                        dbc.Label('Select Secondary Histogram Color:'),
                        dbc.Select(id='secondary_c',
                                   options = [{'label':k, 'value':k} for k in full_color_list],
                                   value='orange'),
                        ])], style={'padding': '10px'}),
                    ]), width=3),
            dbc.Col(dbc.Card([
                dbc.FormGroup([dcc.Loading(id='z-score-graph-loading', children = [html.Div(children = [dcc.Graph(id='output_z_score')], style={"height": "750px", "overflow": "scroll"})], type='default')])
            ]), width=6),
            dbc.Col(                
                    dbc.Card([dbc.FormGroup([html.Div([
                    html.H4("Explore Genes", className="display-5"),
                    html.P('Complete the form on the left to generate a table of predictions of selective genes:', style={'font-style': 'italic'}), 
                    html.Div(id='selective-gene-table',
                             children = [dcc.Loading(id='load_table', children=[dbc.Table.from_dataframe(df_dummy_1[['Gene', 'CERES Pred', 'Z-Score Pred']], striped=True, bordered=True, hover=True, size='sm')], type='default'),], 
                                 style={"height": "580px", "overflow": "scroll",  # enable scrolling
                                        },),
                    html.Br(), 
                    html.Div([dbc.Button('Download Genes', id='download-z-genes-button', color="primary", n_clicks=0, disabled = True,
                                style={'margin-right': '10px'}),
                             dbc.Button('GO', id='run-gprofiler-zscores', color="info", n_clicks=0, disabled = True,
                                style={'margin-left': '10px'}),], 
                    style = {'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),
                    ],
                         )])], style={'padding': '10px'}), width=3),
        ]),

        html.Br(),
        ], style=SIDEBAR_STYLE))]),

volcano_frame = dbc.Row([dbc.Col(html.Div([
        html.P('Select an experiment below to visualize predictions of significant genes based on gene effect score and p-value thresholds', 
        style={'font-size': '1.2rem', 'font-style': 'italic'}),

        dbc.Row([
            dbc.Col(html.Div([
                    dbc.Card([dbc.FormGroup([
                        html.H6('Data:', className="card-title"),
                        dcc.Dropdown(id='choose_experiment_4',
                                   placeholder='select an experiment...'),
                        html.Br(),
                        dcc.Dropdown(id='choose_gene_cats',
                                   placeholder='select gene categories...',
                                   options = [{"label":'Conditional Essential', "value": "conditional essential"},
                                                {"label":'Nonessential', "value": "common nonessential"},
                                                {"label":'Essential', "value": "common essential"},
                                                ],
                                    multi=True),
                        html.Br(),
                        html.H6('Colors:', className="card-title"),
                        dcc.Dropdown(id='conditions_c',
                                     placeholder='select colormap for gene categories...',
                                     options = [{'label':k, 'value':k} for k in colormap_list],
                                     value = 'rainbow'),
                        html.Br(),
                        dcc.Dropdown(id='significant_c',
                                     placeholder='choose color for significant genes...',
                                     options = [{'label':k, 'value':k} for k in full_color_list],
                                     value='red'),
                        html.Br(),
                        dcc.Dropdown(id='gene_effect_c',
                                     placeholder='choose color for gene effect threshold...',
                                     options = [{'label':k, 'value':k} for k in full_color_list],
                                     value='purple'),
                        html.Br(),
                        dcc.Dropdown(id='p_val_c',
                                     placeholder='choose color for p-val threshold...',
                                     options = [{'label':k, 'value':k} for k in full_color_list],
                                     value='maroon'),
                        html.Br(),
                        html.H6('Thresholds:', className="card-title"),
                        dbc.Label('select CERES thresholds:'),
                        dcc.RangeSlider(-2, 2, 0.1, value=[-0.5, 0.5], id='ceres_range_slider',
                        marks=dict(zip([i for i in np.arange(-2.5, 2.5, 0.5)], ["{:.1f}".format(i) for i in np.arange(-2.5, 2.5, 0.5)])),
                        tooltip={"placement": "bottom", "always_visible": False}
                        ),
                        #html.Br(),
                        dbc.Label('select p-val threshold:'),
                        dcc.Slider(0, 10, 1, value=2, id='p_val_threshold'),
                        ])], style={'padding': '10px'}),
                    ]), width=3),
            dbc.Col(dbc.Card([
                dbc.FormGroup([dcc.Loading(id='volcano-graph-loading', children = [html.Div(children = [dcc.Graph(id='output_volcano_plot')], style={"height": "750px", "overflow": "scroll"})], type='default')])
            ]), width=6),
            dbc.Col(                
                    dbc.Card([dbc.FormGroup([html.Div([
                    html.H4("Explore Genes", className="display-5"),
                    html.P('Complete the form on the left to generate a table of predictions of significant selective genes:', style={'font-style': 'italic'}), 
                    html.Div(id='significant-gene-table-p',
                             children = [dcc.Loading(id='load_table', children=[dbc.Table.from_dataframe(df_dummy_1[['Gene', 'CERES Pred', 'Z-Score Pred']], striped=True, bordered=True, hover=True, size='sm')], type='default'),], 
                                 style={"height": "580px", "overflow": "scroll",  # enable scrolling
                                        },),
                    html.Br(), 
                    html.Div([dbc.Button('Download Genes', id='download-p-genes-button', color="primary", n_clicks=0, disabled = True,
                                style={'margin-right': '10px'}),
                             dbc.Button('GO', id='run-gprofiler-pvalue', color="info", n_clicks=0, disabled = True,
                                style={'margin-left': '10px'}),], 
                    style = {'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),
                    ],
                         )])], style={'padding': '10px'}), width=3),
        ]),

        html.Br(),
        ], style=SIDEBAR_STYLE))]),


layout = dbc.Container([
    html.Div([
    html.Br(),
    dcc.Loading(id='warning_overview_load', children = [html.Div([html.Label(id='data_available_overview', style={'color':'#548235'})], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})], type='default'),
    
    dcc.Tabs(id="tabs-styled-with-inline-1", value='distributions-tab', children = [
        dcc.Tab(label='''Overview''', value='distributions-tab', style=tab_style, selected_style=tab_selected_style, children=distributions_frame),
        dcc.Tab(label='''Volcano Plot''', value='volcano-tab', style=tab_style, selected_style=tab_selected_style, children=volcano_frame),
        dcc.Tab(label='''Z-Score Hits''', value='z-score-tab', style=tab_style, selected_style=tab_selected_style, children=z_score_frame),
        dcc.Tab(label='''Table''', value='table-tab', style=tab_style, selected_style=tab_selected_style, children=table_frame),
    ]),

    html.Div(id='dummy_div_overview'),
    html.Div(id='dummy_div_2'),], 

    style={
        "transform": "scale(0.85)",
        "transform-origin": "top",
        "className": 'small'
    },

             ),
], fluid=True, style={'className': 'small'})


@callback(Output('choose_experiment_1', 'options'),
          Output('choose_experiment_2', 'options'),
          Output('choose_experiment_3', 'options'),
          Output('choose_experiment_4', 'options'),
          State('experiment_labels', 'data'),
          Input('dummy_div_overview', 'children'),)
def update_dropdown(experiments, aux):
    if experiments is None:
        raise PreventUpdate
    else:
        return [{'label': e, 'value': e} for e in experiments], [{'label': e, 'value': e} for e in experiments], [{'label': e, 'value': e} for e in experiments], [{'label': e, 'value': e} for e in experiments]

@callback(Output('data_available_overview', 'children'),
          State('df_pred-store', 'data'),
          State('experiment_labels', 'data'),
          Input('dummy_div_overview', 'children'),)
def update_dropdown(experiments, pred_data, aux):
    if experiments is None or pred_data is None:
        return [dbc.Alert(
                    [html.Div([
                    html.I(className="bi bi-info-circle-fill me-2"),
                    html.Br(),
                    dbc.ListGroup(
                    [dbc.ListGroupItem(['Prediction data is unavailable - Please visit the "Predict" tab to run a new inference or load a sample prediction'])],),
                    ],),
                    ], color="danger",
                    className="d-flex align-items-center",
                    style={'align-items': 'center'}
                    ),]
    else:
        return []