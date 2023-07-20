import dash
from dash import html, callback
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import pickle
import pandas as pd

df_clusters = pd.read_csv('./data/feat_summary_varExp_filtered_class.csv')
landmark_genes = sorted(list(set(df_clusters['feature'].tolist())))

dash.register_page(__name__)

SIDEBAR_STYLE = {
    "padding": "1rem 1rem",
    "background-color": "#f8f9fa",
}

df_dummy_1 = pd.DataFrame(columns=['Gene', 'Gene Category', 'Average', 'Standard Deviation', 'CERES Pred', 'Z-Score Pred'])

layout = dbc.Container([

    html.Div([
    html.Br(),

    dbc.Row([dbc.Col(html.Div([
        html.H4('Community Analysis', className="display-5"),
        html.P('Use this tool to visualize relevant functional connections between different genes within a Louvain network for a specific cell-line or experiment.'),

        dbc.Row([dbc.Col(html.Div([
                   html.Br(),
                   dbc.Card([dbc.FormGroup([html.H6('Community:', className="card-title"),
                                            dbc.Label('Select a Louvain Community by Landmark Feature'),
                                            dbc.Select(
                                               id="select_community",
                                               options=[{'label': value, 'value': value} for value in landmark_genes],
                                               #labelStyle={'display': 'inline-block', 'width': '12em', 'line-height':'0.5em'}
                                               ),
                                               html.Br(),
                                               html.Br(),
                                               html.H6('Cell-Line:', className="card-title"),
                                                dbc.Label('Select a Cell-Line:'),
                                                dbc.Select(
                                               id="select_experiment",
                                               #labelStyle={'display': 'inline-block', 'width': '12em', 'line-height':'0.5em'}
                                               ),
                                               html.Br(),
                                               html.Br(),
                                            dbc.Label('or Select an Experiment:'),
                                            dbc.Select(
                                                id="select_experiment_cluster"
                                            ),
                                                html.Br(),
                                               html.Br(),
                                                html.H6('Layout:', className="card-title"),
                                                dbc.Label('Select a Network Layout:'),
                                                dbc.Select(
                                               id="select_network_layout",
                                               options=[{'label':'  Spring', 'value':'spring'},
                                                        {'label':'  Circular', 'value':'circular'},
                                                        {'label':'  Random', 'value':'random'},
                                                        {'label':'  Spectral', 'value':'spectral'},
                                                        {'label':'  Fruchterman Reingold', 'value':'fruchterman_reingold'},],
                                               #labelStyle={'display': 'inline-block', 'width': '12em', 'line-height':'0.5em'}
                                               ),

                                               html.Br(),
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
                 dbc.Col(
            html.Div([
                        html.P(children = '', id='network_graph_1_name', style={'textAlign': 'center'}),
                        dcc.Loading(id='network_1_load', children = [dcc.Graph(id='network_graph_1',  figure = {"layout": {"height": 700}})], type='default'),],
                             ), width=6
            ), 
                dbc.Col(                
                    dbc.Card([dbc.FormGroup([html.Div([
                    html.H4("Explore Genes", className="display-5"),
                    html.P('Complete the form to the left to generate a table of genes within the louvain community selected:'),
                    html.Br(),   
                    html.Div(id='louvain-gene-table',
                             children = [dcc.Loading(id='load_table', children=[dbc.Table.from_dataframe(df_dummy_1[['Gene', 'CERES Pred', 'Z-Score Pred']], striped=True, bordered=True, hover=True)], type='default'),], 
                                 style={
                                            "height": "550px",
                                            "overflow": "scroll",  # enable scrolling
                                        },),
                    ],
                         )])], style={'padding': '10px'}), width=3)
              ]),
        ], style=SIDEBAR_STYLE))]),

        html.Br(),

    dbc.Row([dbc.Col(html.Div([
        html.H4('Cluster Analysis', className="display-5"),
        html.P('Use this tool to visualize PHATE clusters within a Louvain community and find genes which define that cluster.'),

        dbc.Row([dbc.Col(html.Div([
                   html.Br(),
                   dbc.Card([dbc.FormGroup([html.H6('Community:', className="card-title"),
                                            dbc.Label('Select a Louvain Community by Landmark Feature'),
                                            dbc.Select(
                                               id="select_community_2",
                                               options=[{'label': value, 'value': value} for value in landmark_genes],
                                               ),
                                                html.Br(),
                                                html.Br(),
                                                html.H6('Dimensional Reduction:', className="card-title"),
                                                dbc.Label('Select a Visualization:'),
                                                dbc.Select(
                                               id="select_cluster_layout",
                                               options=[{'label':'  PCA', 'value': '0PCA'},
                                                        {'label':'  tSNE', 'value': '1tSNE'},
                                                        {'label':'  UMAP', 'value': '3UMAP'},
                                                        {'label':'  PHATE', 'value': '2PHATE'},
                                                        {'label':'  PHATE (tSNE)', 'value': '4tSNE'},
                                                        {'label':'  PHATE (UMAP)', 'value': '5UMAP'},],
                                               ),
                                                html.Br(),
                                                html.Br(),
                                                html.H6('Cluster:', className="card-title"),
                                                dbc.Label('Select a Cluster:'),
                                                dbc.Select(
                                               id="select_cluster_2",
                                               ),
                                               html.Br(),
                                               html.Br(),
                                               html.H6('Cell-Line:', className="card-title"),
                                                dbc.Label('Select a Cell-Line:'),
                                                dbc.Select(
                                               id="select_cell_cluster_2",
                                               ),
                                               html.Br(),
                                               html.Br(),
                                               dbc.Label('or Select an Experiment:'),
                                               dbc.Select(
                                                id="select_experiment_cluster_2"
                                            ),
                                               ])],
                            style={'padding': '10px'}),
                    
            ]), width=3), 
                 dbc.Col(
            html.Div([
                        html.P(children = '', id='cluster_graph_2_name', style={'textAlign': 'center'}),
                        dcc.Loading(id='network_2_load', children = [dcc.Graph(id='cluster_graph_1',  figure = {"layout": {"height": 700}})], type='default'),],
                             ), width=6
            ), 
                dbc.Col(                
                    dbc.Card([dbc.FormGroup([html.Div([
                    html.H4("Cluster Features", className="display-5"),
                    html.P('Complete the form to the left to generate a table of top features for a selected cluster and cell-line:'),
                    html.Br(),   
                    html.Div(id='cluster-gene-table',
                             children = [dcc.Loading(id='load_table_2', children=[dbc.Table.from_dataframe(df_dummy_1[['Gene', 'CERES Pred', 'Z-Score Pred']], striped=True, bordered=True, hover=True)], type='default'),], 
                                 style={
                                            "height": "550px",
                                            "overflow": "scroll",  # enable scrolling
                                        },),
                    ],
                         )])], style={'padding': '10px'}), width=3)
              ]),

        html.Br(),

        dbc.Row([dbc.Col(html.Div([
            html.H4('Cluster Feature Importance', className="display-5"),
            html.P('Select Clusters to Visualize Most Important Features:')]))]),

        dbc.Row([dbc.Col(html.Div([
                   html.Br(),
                    dbc.Select(
                    id="select_clust_feats_1",),
                    html.Br(),
                    html.Br(),
                    dcc.Loading(id='clust_feats_1_load', children = [dcc.Graph(id='cluster_feats_1',  figure = {"layout": {"height": 700}})], type='default'),
            ]), width=3), 
                dbc.Col(html.Div([
                   html.Br(),
                    dbc.Select(
                    id="select_clust_feats_2",),
                    html.Br(),
                    html.Br(),
                    dcc.Loading(id='clust_feats_2_load', children = [dcc.Graph(id='cluster_feats_2',  figure = {"layout": {"height": 700}})], type='default'),
            ]), width=3), 
                dbc.Col(html.Div([
                   html.Br(),
                    dbc.Select(
                    id="select_clust_feats_3",),
                    html.Br(),
                   html.Br(),
                    dcc.Loading(id='clust_feats_3_load', children = [dcc.Graph(id='cluster_feats_3',  figure = {"layout": {"height": 700}})], type='default'),
            ]), width=3), 
                dbc.Col(html.Div([
                   html.Br(),
                    dbc.Select(
                    id="select_clust_feats_4",),
                    html.Br(),
                   html.Br(),
                    dcc.Loading(id='clust_feats_4_load', children = [dcc.Graph(id='cluster_feats_4',  figure = {"layout": {"height": 700}})], type='default'),
            ]), width=3), 
              ]),
        ], style=SIDEBAR_STYLE))]),

        html.Div(id='dummy_div_genes'),],
    
        style={
        "transform": "scale(0.85)",
        "transform-origin": "top",
    },
    )
], fluid = True 
+1, style={'className': 'small'})