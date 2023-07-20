import dash
from dash import html, callback
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
#from dash_extensions.enrich import Dash, DashProxy, Output, Input, State, ServersideOutput, html, dcc, ServersideOutputTransform
import pandas as pd

dash.register_page(__name__, path='/')

SIDEBAR_STYLE = {
    "padding": "1rem 1rem",
    "background-color": "#f8f9fa",
}

layout = dbc.Container([
    html.Br(),
    
    dbc.Row([dbc.Col(html.Div([
        html.H4('GO-LoCo', className="display-4", style={'textAlign': 'center'}),
        html.P('A webapp to perform genome wide loss of function inferences with compressed subsets of genes powered by Lossy Compression (LoCo)', style={'textAlign': 'center'}),
        html.Br(),
        ],style=SIDEBAR_STYLE))]),

    html.Br(),

    html.Div(id='dummy_div_homepage'),
])
