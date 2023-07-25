import dash
from dash import html, callback
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import pickle
import pandas as pd

DepMap19Q4_GENES_PICKLE = './data/19q4_genes.pkl'

dash.register_page(__name__)

SIDEBAR_STYLE = {
    "padding": "1rem 1rem",
    "background-color": "#f8f9fa",
}

layout = dbc.Container([

    html.Div([
    html.Br(),
    
    dbc.Row([dbc.Col(html.Div([
        html.H4('Single Gene Comparisons', className="display-5"),
        html.P('Select a gene below for each plot to visualize the predicted CERES score for all experiments in the backdrop of the distribution of CERES scores for that gene in all other cell lines'),
        dbc.Row([dbc.Col(
                        dcc.Loading(id='louv_load_1', children=[dcc.Dropdown(id='pick-community-1', searchable = True, placeholder='select a louvain community')], type='default'), 
                        width=6)]),
        html.Br(),
        dbc.Row([dbc.Col(
            dcc.Loading(id='gene_load', children=[dcc.Dropdown(id='gene-list-dropdown_1', searchable = True, placeholder='select a gene')], type='default')
            ),
                 dbc.Col(
                     dcc.Loading(id='gene_load', children=[dcc.Dropdown(id='gene-list-dropdown_2', searchable = True, placeholder='select a gene')], type='default')
                     ),
                 ]),
        html.Br(),
        dbc.Row([dbc.Col(
            html.Div([
                        html.P(children = '', id='single_dist_graph_1_name', style={'textAlign': 'center'}),
                        dcc.Graph(id='single_dist_graph_1'),],
                             ), width=6
            ),
            
                 dbc.Col(
                     html.Div([
                        html.P(children = '', id='single_dist_graph_2_name', style={'textAlign': 'center'}),
                        dcc.Graph(id='single_dist_graph_2'),],
                             ), width=6
                     ),]),
        html.Br(),
        dbc.Row([dbc.Col(
            html.Div([
                dbc.Row([dbc.Col(width=1), dbc.Col(
                    html.Div(children = [], id = 'depmap_link_1', style={'textAlign': 'center'})
                
                ),
                        dbc.Col(
                            html.Div(children = [], id = 'uniprot_link_1', style={'textAlign': 'center'})

                        ),
                        dbc.Col(
                            html.Div(children = [], id = 'tumor_portal_link_1', style={'textAlign': 'center'})

                        ),
                        dbc.Col(
                            html.Div(children = [], id = 'ncbi_link_1', style={'textAlign': 'center'})

                        ), dbc.Col(width=1),
                ])
                

            ])
        ),
                dbc.Col(
                    html.Div([
                    dbc.Row([dbc.Col(width=1), dbc.Col(
                        html.Div(children = [], id = 'depmap_link_2', style={'textAlign': 'center'})
                
                        ),
                        dbc.Col(
                            html.Div(children = [], id = 'uniprot_link_2', style={'textAlign': 'center'})

                        ),
                        dbc.Col(
                            html.Div(children = [], id = 'tumor_portal_link_2', style={'textAlign': 'center'})

                        ),
                        dbc.Col(
                            html.Div(children = [], id = 'ncbi_link_2', style={'textAlign': 'center'})

                        ), dbc.Col(width=1),
                ])

                    ])
                )
        ])

        ], style=SIDEBAR_STYLE)
        )]),

    html.Br(),

    dbc.Row([dbc.Col(html.Div([
        html.H4('Single Pairwise Genes Regression Analysis', className="display-5"),
        html.P('Select two genes below to visualize the predicted CERES scores against one another and amongst all other cell lines'),
                
        dbc.Row([dbc.Col(
            dcc.Loading(id='louv_load_1', children=[dcc.Dropdown(id='pick-community-2', searchable = True, placeholder='select a louvain community')], type='default'), 
        width=6)]),

        html.Br(),
        dbc.Row([dbc.Col(
            dcc.Loading(id='gene_load', children=[dcc.Dropdown(id='gene-list-dropdown_3', searchable = True, placeholder='select first gene')], type='default')
            ),
                 dbc.Col(
                     dcc.Loading(id='gene_load', children=[dcc.Dropdown(id='gene-list-dropdown_4', searchable = True, placeholder='select second gene')], type='default')
                     ),
                 ]),
        html.Br(),
        dbc.Row([dbc.Col(html.Div([
                   html.Br(),
                   dbc.Card([dbc.FormGroup([html.H6('Color:', className="card-title"),
                                            dbc.Label('Select Color By:'),
                                            dbc.Select(
                                               id="select_color_by",
                                               options=[{'label':'  Lineage', 'value':'lineage'},
                                                        {'label':'  Lineage Subtype', 'value':'lineage_subtype'},
                                                        {'label':'  Collection Site', 'value':'sample_collection_site'},
                                                        {'label':'  Disease', 'value':'disease'},
                                                        {'label':'  Disease Subtype', 'value':'disease_subtype'},
                                                        {'label':'  Origin', 'value':'primary_or_metastasis'},],
                                               #labelStyle={'display': 'inline-block', 'width': '12em', 'line-height':'0.5em'}
                                               ),])],
                            style={'padding': '10px'}),
                   html.Br(),
                   dbc.Card([dbc.FormGroup([html.H6('Linear Regression:', className="card-title"),
                                            dbc.Label('Select Subcategory:'),
                                            dbc.Select(
                                               id="select_category",
                                               #labelStyle={'display': 'inline-block', 'width': '12em', 'line-height':'0.5em'}
                                               ),
                                            ])],
                            style={'padding': '10px'}),                   
            ]), width=3), 
                 dbc.Col(
            html.Div([
                        html.P(children = '', id='pairwise_gene_graph_name', style={'textAlign': 'center'}),
                        dcc.Loading(id='pairwise_load', children = [dcc.Graph(id='pairwise_gene_comp',  figure = {"layout": {"height": 700}})], type='default'),],
                             ), width=9
            )
                 ]),
        ], style=SIDEBAR_STYLE))]),

    html.Br(),

    dbc.Row(
        dbc.Col(
                html.Div([
                    html.H4('Multiple Pairwise Genes Regression Analysis', className="display-5"),
                    html.P('Select between 2 to 10 genes in the dropdown below to visualize the pairwise comparisons of all genes across all cell lines including your predicted CERES scores'),
                    
                    dbc.Row([dbc.Col(
                        dcc.Loading(id='louv_load_1', children=[dcc.Dropdown(id='pick-community-3', searchable = True, placeholder='select a louvain community')], type='default'), 
                    width=6)]),

                    html.Br(),
                    dbc.Row(dbc.Col(dcc.Dropdown(id='multi-gene-list-dropdown', searchable = True, placeholder='select second gene', multi = True))),
                    html.Br(),
                    dbc.Row([dbc.Col(), dbc.Col(html.Div([
                        html.P(children = '', id='multi_gene_graph_name', style={'textAlign': 'center'}),
                        dcc.Loading(id='multi-gene-comp', children=[dcc.Graph(id='multi_gene_comp', figure = {"layout": {"height": 1000}})], type='default'),],
                             ), width = 10), dbc.Col()])], style=SIDEBAR_STYLE)
                )),

    html.Div(id='dummy_div_genes'),],

    style={
        "transform": "scale(1)",
        "transform-origin": "top",
    },
             ),

], fluid = True, style={'className': 'small'})

'''@callback(Output('gene-list-dropdown_1', 'options'),
          Output('gene-list-dropdown_2', 'options'),
          Output('gene-list-dropdown_3', 'options'),
          Output('gene-list-dropdown_4', 'options'),
          Output('pick-community-1', 'options'),
          Output('pick-community-2', 'options'),
          Output('pick-community-3', 'options'),
          Input('dummy_div_genes', 'children'),)
def update_dropdown(aux):
    with open(DepMap19Q4_GENES_PICKLE, 'rb') as f:
        gene_list = pickle.load(f)

    df_clusters = pd.read_csv('./data/feat_summary_varExp_filtered_class.csv')
    landmark_genes = sorted(list(set(df_clusters['feature'].tolist())))
    
    return gene_list, gene_list, gene_list, gene_list, landmark_genes, landmark_genes, landmark_genes'''