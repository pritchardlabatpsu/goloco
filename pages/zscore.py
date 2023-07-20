import dash
from dash import html, callback
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import pandas as pd
#from dash_extensions.enrich import Dash, DashProxy, Output, Input, State, ServersideOutput, html, dcc, ServersideOutputTransform

dash.register_page(__name__)

SIDEBAR_STYLE = {
    "padding": "1rem 1rem",
    "background-color": "#f8f9fa",
}

layout = dbc.Container([
    html.Br(),
    
    dbc.Row([dbc.Col(html.Div([
        html.H4('Predicted Z-Scores', className="display-5"),
        html.P('Select an experiment below and select gene categories to visualize CERES Score proedictions by Z-Scores'),
        html.Div([dcc.Dropdown(id='choose_experiment_3')], style={"width": "40%"}),
        html.Br(),
        dcc.Checklist([{"label":'Conditional Essential', "value": "conditional essential"},
                       {"label":'Nonessential', "value": "common nonessential"},
                       {"label":'Essential', "value": "common essential"},
                       ], inline=True, inputStyle={"margin-right": "6px", "margin-left": "15px"},
                      id='gene_cat_checklist_2'),
        html.Br(),
        html.P('CERES Predictions by Z-Scores', style={'textAlign': 'center'}),
        dcc.Loading(id='z-score-graph-loading', children = [html.Div([dcc.Graph(id='output_z_score')])], type='default'),
        ],style=SIDEBAR_STYLE))]),

    html.Br(),

    html.Div(id='dummy_div_zscore'),
], fluid = True, style={'padding-right': '100px',
                      'padding-left': '100px',
                      'padding-top': '20px'})

@callback(Output('choose_experiment_3', 'options'),
          State('experiment_labels', 'data'),
          Input('dummy_div_zscore', 'children'),)
def update_dropdown(experiments, aux):
    if experiments is None:
        raise PreventUpdate
    else:
        return experiments