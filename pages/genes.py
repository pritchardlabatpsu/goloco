import dash
from dash import html, callback
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import pickle
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import matplotlib

DepMap19Q4_GENES_PICKLE = './data/19q4_genes.pkl'

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
    'backgroundColor': '#0f8de5',  #'#119DFF'
    'color': 'white',
    'font-family': 'Helvetica',
    'fontSize': '20px',
    'padding': '6px',
    'width': '20%',
    'borderWidth': '1px',
    'borderStyle': 'solid',
    'borderRadius': '1px',
}

single_gene_comp = html.Div([
        html.P('Select genes below to visualize the predicted CERES scores across experiments compared to the distribution of scores across all other DepMap cell lines',
        style={'font-size': '1.2rem', 'font-style': 'italic'}),
        dbc.Row([dbc.Col(
                        dcc.Loading(id='louv_load_1', children=[dcc.Dropdown(id='pick-community-1', searchable = True, placeholder='Select a Louvain community to narrow gene selection...')], type='default'), 
                        width=6)]),
        html.Br(),

        dbc.Card([dbc.FormGroup([
        dbc.Row([dbc.Col(
            html.Div([  dcc.Loading(id='gene_load', children=[dcc.Dropdown(id='gene-list-dropdown_1', searchable = True, placeholder='Select a gene...')], type='default'),
                        html.P(children = '', id='single_dist_graph_1_name', style={'textAlign': 'center'}),
                        dcc.Graph(id='single_dist_graph_1'),],
                             ), width=6
            ),
            
                 dbc.Col(
                     html.Div([
                        dcc.Loading(id='gene_load', children=[dcc.Dropdown(id='gene-list-dropdown_2', searchable = True, placeholder='Select a gene...')], type='default'),
                        html.P(children = '', id='single_dist_graph_2_name', style={'textAlign': 'center'}),
                        dcc.Graph(id='single_dist_graph_2'),],
                             ), width=6
                     ),]),
        html.Br(),
        dbc.Row([dbc.Col(
            html.Div([
                dbc.Row([dbc.Col(width=1), 
                        dbc.Col(html.Div(children = [], id = 'depmap_link_1', style={'textAlign': 'center'})),
                        dbc.Col(html.Div(children = [], id = 'uniprot_link_1', style={'textAlign': 'center'})),
                        dbc.Col(html.Div(children = [], id = 'tumor_portal_link_1', style={'textAlign': 'center'})),
                        dbc.Col(html.Div(children = [], id = 'ncbi_link_1', style={'textAlign': 'center'})), 
                        dbc.Col(width=1),
                        ])
                ])),

                dbc.Col(
                    html.Div([
                    dbc.Row([dbc.Col(width=1), 
                            dbc.Col(html.Div(children = [], id = 'depmap_link_2', style={'textAlign': 'center'})),
                            dbc.Col(html.Div(children = [], id = 'uniprot_link_2', style={'textAlign': 'center'})),
                            dbc.Col(html.Div(children = [], id = 'tumor_portal_link_2', style={'textAlign': 'center'})),
                            dbc.Col(html.Div(children = [], id = 'ncbi_link_2', style={'textAlign': 'center'})), 
                            dbc.Col(width=1),])
                    ]))])
        ])], style={'padding': '10px'}),

                    ], style=SIDEBAR_STYLE)

pairwise_gene_comp = html.Div([
        html.P('Select two genes below to visualize the predicted CERES Scores against one another and among all other DepMap cell lines',
        style={'font-size': '1.2rem', 'font-style': 'italic'}),
        dbc.Row([dbc.Col(html.Div([
                    dbc.Card([dbc.FormGroup([
                    html.H6('Select Genes:', className="card-title"),
                    dbc.Label('Louvain Community:'),
                    dcc.Loading(id='louv_load_1', children=[dcc.Dropdown(id='pick-community-2', searchable = True)], type='default'), 
                    html.Br(),
                    dbc.Label('Select First Gene:'),
                    dcc.Loading(id='gene_load', children=[dcc.Dropdown(id='gene-list-dropdown_3', searchable = True)], type='default'),
                    html.Br(),
                    dbc.Label('Select Second Gene:'),
                    dcc.Loading(id='gene_load', children=[dcc.Dropdown(id='gene-list-dropdown_4', searchable = True)], type='default'),
                    html.Br(),
                    html.H6('Color:', className="card-title"),
                    dbc.Label('Select Color By:'),
                            dbc.Select(
                                id="select_color_by",
                                options=[{'label':'  Lineage', 'value':'lineage'},
                                        {'label':'  Lineage Subtype', 'value':'lineage_subtype'},
                                        {'label':'  Collection Site', 'value':'sample_collection_site'},
                                        {'label':'  Disease', 'value':'disease'},
                                        {'label':'  Disease Subtype', 'value':'disease_subtype'},
                                        {'label':'  Origin', 'value':'primary_or_metastasis'},],
                                ),
                   html.Br(),
                   html.Br(),
                   html.H6('Linear Regression:', className="card-title"),
                   dbc.Label('Select Subcategory:'),
                   dbc.Select(id="select_category"),]
                   )], style={'padding': '10px'}),
            ]), width=3), 
                 dbc.Col(
                    dbc.Card([dbc.FormGroup([
            html.Div([
                        html.P(children = '', id='pairwise_gene_graph_name', style={'textAlign': 'center'}),
                        dcc.Loading(id='pairwise_load', children = [dcc.Graph(id='pairwise_gene_comp',  figure = {"layout": {"height": 600}})], type='default'),],
                             )])], style={'padding': '10px'}), width=9
            )
                 ]),
        ], style=SIDEBAR_STYLE)

multi_gene_comp = html.Div([
                html.P('Select between 2 and 6 genes in the dropdown below to visualize the pairwise comparisons of all genes across predictions and DepMap cell lines',
                style={'font-size': '1.2rem', 'font-style': 'italic'}),
                html.Br(),
                dbc.Card([dbc.FormGroup([
                    html.H6('Select Genes:', className='card-title'),
                    dbc.Row([
                        dbc.Col(dcc.Loading(id='louv_load_1', children=[dcc.Dropdown(id='pick-community-3', searchable = True, placeholder='Select a Louvain community to narrow gene selection...')], type='default'), ),
                        dbc.Col(),
                    ]),
                    html.Br(),
                    dcc.Dropdown(id='multi-gene-list-dropdown', searchable = True, placeholder='Select genes...', multi = True),
                    html.Br(),
                    html.H6('Primary Scatter Options:', className='card-title'),
                    dbc.Row([
                        dbc.Col(dcc.Dropdown(id='m_primary_color', placeholder='Select a color for secondary scatter...', options = [{'label':k, 'value':k} for k in full_color_list], value='magenta')),
                        dbc.Col(dcc.Dropdown(id='m_experiment_color', placeholder='Select a colormap...', options = [{'label':k, 'value':k} for k in colormap_list], value='viridis')),
                        dbc.Col(),
                    ]),
                    html.Br(),
                    html.H6('Secondary Scatter Options:', className='card-title'),
                    dbc.Row([
                        dbc.Col(dcc.Dropdown(id="m_select_secondary_cat",
                                    placeholder='Select a category...',
                                    options=[{'label':'  Lineage', 'value':'lineage'},
                                            {'label':'  Lineage Subtype', 'value':'lineage_subtype'},
                                            {'label':'  Collection Site', 'value':'sample_collection_site'},
                                            {'label':'  Disease', 'value':'disease'},
                                            {'label':'  Disease Subtype', 'value':'disease_subtype'},
                                            {'label':'  Origin', 'value':'primary_or_metastasis'},],)),
                        dbc.Col(dcc.Dropdown(id="m_select_sub_cat", placeholder='Select a category to load options...'),),
                        dbc.Col(dcc.Dropdown(id='m_secondary_color', placeholder='Select a color for secondary scatter...', options = [{'label':k, 'value':k} for k in full_color_list], value='forestgreen')),
                    ],),
                    html.P(children = '', id='multi_gene_graph_name', style={'textAlign': 'center'}),
                    dbc.Row([
                        dbc.Col(dcc.Loading(id='multi-gene-comp', children=[dcc.Graph(id='multi_gene_comp', figure = {"layout": {"height": 1000}})], type='default'), width=8),
                        dbc.Col(dcc.Loading(id='multi-gene-comp-2', children=html.Div(id='multi_gene_comp_network', children=[]),),),
                    ]),
                ])], style={'padding': '10px'}),
                        
                ], style=SIDEBAR_STYLE)

layout = dbc.Container([
    html.Div([
    html.Br(),
    dcc.Loading(id='warning_genes_load', children = [html.Div([html.Label(id='data_available_genes', style={'color':'#548235'})], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})], type='default'),

    dcc.Tabs(id="tabs-styled-with-inline-1", value='multi-lin-reg-tab', children = [
        dcc.Tab(label='''Multi Linear Regression''', value='multi-lin-reg-tab', style=tab_style, selected_style=tab_selected_style, children=multi_gene_comp),
        dcc.Tab(label='''Pairwise Linear Regression''', value='linear-reg-tab', style=tab_style, selected_style=tab_selected_style, children=pairwise_gene_comp),
        dcc.Tab(label='''Gene Distributions''', value='gene-dist-tab', style=tab_style, selected_style=tab_selected_style, children=single_gene_comp),
    ]),

    html.Div(id='dummy_div_genes'),],

    style={
        "transform": "scale(0.85)",
        "transform-origin": "top",
        "className": 'small'
    },
             ),

], fluid = True, style={'className': 'small'})

@callback(Output('data_available_genes', 'children'),
          State('df_pred-store', 'data'),
          State('experiment_labels', 'data'),
          Input('dummy_div_genes', 'children'),)
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