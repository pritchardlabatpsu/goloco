import dash
from dash import html, callback
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import pickle
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib


cell_line_DIRECTORY = './data/cell_line_names.pkl'

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

df_dummy_1 = pd.DataFrame(columns=['Gene', 'Gene Category', 'Average', 'Standard Deviation', 'CERES Pred', 'Z-Score Pred'])

community_layout = html.Div([
        html.P('Use this tool to visualize functional connections between different genes within a Louvain community for a DepMap cell line or experiment', 
        style={'font-size': '1.2rem', 'font-style': 'italic'}),

        dbc.Row([dbc.Col(html.Div([
                   dbc.Card([dbc.FormGroup([html.H6('Network:', className="card-title"),
                                            dcc.Dropdown(id="select_community",
                                               options=[{'label': 'louvain_' + str(i), 'value': i} for i in range(44)],
                                               multi=True,
                                               searchable=True,
                                               placeholder='select a louvain community'),
                                                html.Br(),
                                                dcc.Dropdown(id="select_network_layout",
                                               options=[{'label':'  Spring', 'value':'cose'},
                                                        {'label':'  Circular', 'value':'circle'},
                                                        {'label':'  Grid', 'value':'grid'},
                                                        {'label':'  Random', 'value':'random'},
                                                        {'label':'  Concentric', 'value':'concentric'},
                                                        {'label':'  Breadthfirst', 'value':'breadthfirst'},],
                                                placeholder='select a network layout',
                                                value='cose',
                                               ),
                                            html.Br(),
                                            html.H6('Data:', className="card-title"),
                                            dcc.Dropdown(id="select_experiment", placeholder='select a cell line...'),
                                            html.Br(),
                                            dcc.Dropdown(id="select_experiment_cluster", placeholder='...or select an experiment'),
                                            html.Br(),
                                            html.H6('Nodes:', className="card-title"),
                                            dcc.Dropdown(id="select_node_color_by", placeholder='select node color by..',
                                                options=[{'label':'  Louvain Community', 'value':'community'},
                                                        {'label':'  CERES Score', 'value':'ceres_score'},
                                                        {'label':'  Z Score', 'value':'z_score'},
                                                        {'label':'  Betweeness Centrality', 'value':'betweenness_centrality'},
                                                        {'label':'  Degree Centrality', 'value':'degree_centrality'},
                                                        {'label':'  Closeness Centrality', 'value':'closeness_centrality'},],
                                                value='community'),
                                            html.Br(),
                                            dcc.Dropdown(id='select_node_colormap', placeholder='select a colormap...',
                                                options = [{'label':k, 'value':k} for k in colormap_list],
                                                value='rainbow'),
                                            html.Br(),
                                            dbc.Label('Select Color Resolution:'),
                                            dcc.Slider(1, 10, 1, value=10, id='cres_threshold'),
                                            html.Br(),
                                            dcc.Dropdown(id='select_node_size_by', placeholder='select node size by...',
                                                    options=[{'label':'  CERES Score', 'value':'ceres_score'},
                                                        {'label':'  Z Score', 'value':'z_score'},
                                                        {'label':'  Betweeness Centrality', 'value':'betweenness_centrality'},
                                                        {'label':'  Degree Centrality', 'value':'degree_centrality'},
                                                        {'label':'  Closeness Centrality', 'value':'closeness_centrality'},],
                                                        value='degree_centrality',),
                                            html.Br(),
                                            html.H6('Gene Effect Threshold:', className="card-title"),
                                            dbc.Label('Select a CERES Threshold:'),
                                            dcc.Slider(-1, 0, 0.1, value=-0.5, id='sig_genes_threshold'),
                                            html.Br(),
                                            dbc.Label('Select two-tail Z-Score Threshold:'),
                                            dcc.Slider(0, 2.5, 0.1, value=0, id='z-score_threshold', 
                                                marks={
                                                    0: '0',
                                                    0.5: '0.5',
                                                    1: '1',
                                                    1.5: '1.5',
                                                    2: '2',
                                                    2.5: '2.5'
                                                }, 
                                                tooltip={"placement": "bottom", "always_visible": False})
                                               ])],
                            style={'padding': '10px'}),
                    
            ]), width=3), 
        dbc.Col(dbc.Card([dbc.FormGroup([html.Div([
                        html.P(children = '', id='network_graph_1_name', style={'textAlign': 'center'}),
                        dcc.Loading(id='network_1_load', children = [html.Div(id='network_graph_1', children=[], style={'height': '700px'})], type='default'),],
                        ),])], style={'padding': '10px'}), width=6), 
        dbc.Col(                
            dbc.Card([dbc.FormGroup([html.Div([
            html.H4("Explore Genes", className="display-5"),
            html.P('Complete the form to the left to generate a table of genes within the louvain community selected:',
            style={'font-style': 'italic'}),
            html.Div(id='louvain-gene-table',
                        children = [dcc.Loading(id='load_table', children=[dbc.Table.from_dataframe(df_dummy_1[['Gene', 'CERES Pred', 'Z-Score Pred']], striped=True, bordered=True, hover=True, size='sm')], type='default'),], 
                            style={
                                    "height": "550px",
                                    "overflow": "scroll",  # enable scrolling
                                },),
            ],
                    )])], style={'padding': '10px'}), width=3)
              ]),
        ], style=SIDEBAR_STYLE)

cluster_layout = html.Div([
        html.P('Use this tool to visualize PHATE clusters within a Louvain community and find genes which define that cluster',
        style={'font-size': '1.2rem', 'font-style': 'italic'}),

        dbc.Row([dbc.Col(html.Div([
                   dbc.Card([dbc.FormGroup([html.H6('Community:', className="card-title"),
                                            dcc.Dropdown(
                                               id="select_community_2",
                                               options=[{'label': 'louvain_' + str(i), 'value': i} for i in range(44)],
                                               placeholder='select a louvain community...',
                                               ),
                                                html.Br(),
                                            dcc.Dropdown(id='select_node_colormap_2', placeholder='select a colormap...',
                                                options = [{'label':k, 'value':k} for k in colormap_list],
                                                value='cool'),
                                                html.Br(),
                                                html.H6('Dimensional Reduction:', className="card-title"),
                                                dcc.Dropdown(
                                                id="select_cluster_layout",
                                               options=[{'label':'  PCA', 'value': '0PCA'},
                                                        {'label':'  tSNE', 'value': '1tSNE'},
                                                        {'label':'  UMAP', 'value': '3UMAP'},
                                                        {'label':'  PHATE', 'value': '2PHATE'},
                                                        {'label':'  PHATE (tSNE)', 'value': '4tSNE'},
                                                        {'label':'  PHATE (UMAP)', 'value': '5UMAP'},],
                                                placeholder='select a visualization...',
                                                value='1tSNE',
                                               ),
                                               html.Br(),
                                               html.H6('Cluster:', className="card-title"),
                                               dcc.Dropdown(
                                               id="select_cluster_2",
                                               placeholder='select a cluster...',
                                               ),
                                               html.Br(),
                                               html.H6('Data:', className="card-title"),
                                               dcc.Dropdown(
                                               id="select_cell_cluster_2",
                                               placeholder='select a DepMap cell line...'
                                               ),
                                               html.Br(),
                                               dcc.Dropdown(
                                                id="select_experiment_cluster_2",
                                                placeholder='...or select an experiment',
                                                ),
                                               ])],
                            style={'padding': '10px'}),
                    
            ]), width=3), 
                 dbc.Col(
                    html.Div([dbc.Card([dbc.FormGroup([
                            html.P(children = '', id='cluster_graph_2_name', style={'textAlign': 'center'}),
                            dcc.Loading(id='network_2_load', children = [dcc.Graph(id='cluster_graph_1',  figure = {"layout": {"height": 670}})], type='default'),],
                             ),], style={'padding': '10px'})]), width=6), 
                dbc.Col(                
                    dbc.Card([dbc.FormGroup([html.Div([
                    html.H4("Cluster Features", className="display-5"),
                    html.P('Complete the form on the left to generate a table of top features for specified cluster and cell line:',
                    style={'font-style': 'italic'}),
                    html.Div(id='cluster-gene-table',
                             children = [dcc.Loading(id='load_table_2', children=[dbc.Table.from_dataframe(df_dummy_1[['Gene', 'CERES Pred', 'Z-Score Pred']], striped=True, bordered=True, hover=True, size='sm')], type='default'),], 
                                 style={
                                            "height": "550px",
                                            "overflow": "scroll",  # enable scrolling
                                        },),
                    ],
                         )])], style={'padding': '10px'}), width=3)
              ]),

        html.Br(),

        dbc.Card([
            dbc.FormGroup([
                    dbc.Row([dbc.Col(html.Div([
                    html.Br(),
                    html.Br(),
                    dcc.Loading(id='clust_feats_1_load', children = [dcc.Graph(id='cluster_feats_1',  figure = {"layout": {"height": 700}})], type='default'),
            ]),),]),
            ])
        ])

        ], style=SIDEBAR_STYLE)

layout = dbc.Container([

    html.Div([
    html.Br(),

    dcc.Tabs(id="tabs-styled-with-inline-2", value='community-tab', children = [
        dcc.Tab(label='''Community Analysis''', value='community-tab', style=tab_style, selected_style=tab_selected_style, children=community_layout),
        dcc.Tab(label='''Cluster Analysis''', value='cluster-tab', style=tab_style, selected_style=tab_selected_style, children=cluster_layout),
    ]),

    html.Div(id='dummy_div_clusters'),],
    
    style={
        "transform": "scale(0.85)",
        "transform-origin": "top",
    },)
], fluid = True, style={'className': 'small'})

@callback(Output('select_cell_cluster_2', 'options'),
          Output('select_experiment', 'options'),
          Input('dummy_div_clusters', 'children'),)
def update_cell_line_dropdown(aux):
    cell_lines = pd.read_pickle(cell_line_DIRECTORY)
    return cell_lines, cell_lines

@callback(Output('select_experiment_cluster', 'options'),
          Output('select_experiment_cluster_2', 'options'),
          State('experiment_labels', 'data'),
          Input('dummy_div_clusters', 'children'),)
def update_dropdown(experiments, aux):
    if experiments is None:
        raise PreventUpdate
    else:
        return [{'label': e, 'value': e} for e in experiments], [{'label': e, 'value': e} for e in experiments]