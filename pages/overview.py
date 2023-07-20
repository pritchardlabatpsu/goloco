import dash
from dash import html, callback
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash import dash_table
from dash.dash_table.Format import Format
import pandas as pd

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
        html.H4('Table of Predictions', className="display-5"),
        html.P('Select an experiment below to generate a table of predictions'),
        html.Div([dcc.Dropdown(id='choose_experiment_1', multi=True)], style={"width": "40%"}),
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
        ],style=SIDEBAR_STYLE))]),

    html.Br(),

    dbc.Row([dbc.Col(html.Div([
        dbc.Row([
            dbc.Col(
                html.Div([
                    html.H4("Distributions of Predictions", className="display-5"),
                    html.P('Select an experiment below to visualize the distribution of CERES scores and Z-Scores'),
                    html.Div([dcc.Dropdown(id='choose_experiment_2')], style={"width": "40%"}),
                    html.Br(),
                    dcc.Checklist([{"label":'Conditional Essential', "value": "conditional essential"},
                       {"label":'Nonessential', "value": "common nonessential"},
                       {"label":'Essential', "value": "common essential"},
                       ], inline=True, inputStyle={"margin-right": "6px", "margin-left": "15px"},
                      id='gene_cat_checklist_1'),
                    ],
                         )),
            ]),

        html.Br(),

        dbc.Row([
            dbc.Col(
                    html.Div([
                        html.P('Distribution of CERES Scores', style={'textAlign': 'center'}),
                        dcc.Loading(id='ceres_distribution', children =[dcc.Graph(id='total_dist_graph')], type='default'),],
                             ), width=6
                ),
            dbc.Col(html.Div([
                    html.P('Distribution of Z-Scores', style={'textAlign': 'center'}),
                    dcc.Loading(id='ceres_distribution', children=[dcc.Graph(id='z_score_dist_graph')], type='default'),],
                             ), width=6
                    ),

            ]),

        html.Br(),

        dbc.Row([
            dbc.Col(
                    html.Div([
                        html.P('PDF of CERES Scores', style={'textAlign': 'center'}),
                        dcc.Loading(id='pdf_distribution', children=[dcc.Graph(id='total_dist_graph_pdf')], type='default'),],
                             ), width=6
                ),
            dbc.Col(html.Div([
                    html.P('PDF of Z-Scores', style={'textAlign': 'center'}),
                    dcc.Loading(id='pdf_distribution_zscore', children=[dcc.Graph(id='z_score_dist_graph_pdf')], type='default'),],
                             ), width=6
                    ),

            ]),

        ],style=SIDEBAR_STYLE), width = 9),

             dbc.Col(html.Div([
        dbc.Row([
            dbc.Col(
                html.Div([
                    html.H4("Explore Genes", className="display-5"),
                    html.P('Click on a bin from one of the histograms on the left to view all genes within that bin'),
                    html.Br(),   
                    html.Div(id='hist-gene-table',
                             children = [dcc.Loading(id='load_table', children=[dbc.Table.from_dataframe(df_dummy_1[['Gene', 'CERES Pred', 'Z-Score Pred']], striped=True, bordered=True, hover=True)], type='default'),], 
                                 style={
                                            "height": "1090px",  # set the height of the div
                                            "overflow": "scroll",  # enable scrolling
                                        },),
                    ],
                         )),
            ]),
        ],style=SIDEBAR_STYLE), width = 3),

    ]),
    html.Div(id='dummy_div_overview'),
    html.Div(id='dummy_div_2'),], 

    style={
        "transform": "scale(0.85)",
        "transform-origin": "top",
    },

             ),
], fluid=True, style={'className': 'small'})


@callback(Output('choose_experiment_1', 'options'),
          Output('choose_experiment_2', 'options'),
          State('experiment_labels', 'data'),
          Input('dummy_div_overview', 'children'),)
def update_dropdown(experiments, aux):
    if experiments is None:
        raise PreventUpdate
    else:
        return experiments, experiments