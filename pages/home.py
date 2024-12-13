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
    'textAlign': 'center',
}

layout = dbc.Container([
    html.Br(),
    
    dbc.Row([dbc.Col(html.Div(id='goloco_home',
        children = [html.H4('goloco', className="display-4", style={'textAlign': 'center'}),
        html.P('A webapp to perform genome wide loss of function inferences with compressed subsets of genes powered by Lossy Compression (LoCo)', style={'textAlign': 'center', 'font-size': '125%', 'font-family': 'Arial'}),
        html.Div([html.A('Download the User Manual to Get Started', id='download_manual_link', href="#")], style = {'textAlign' : 'center'}),
        dcc.Download(id="download_manual_data"),
        html.Br(),
        html.Img(src='./assets/webapp_io.jpg', style={'height':'80%', 'width':'80%'}),],
        style=SIDEBAR_STYLE))]),        

    html.Div(id='dummy_div_homepage'),
])

@callback(
    Output("download_manual_data", "data"),
    Input("download_manual_link", "n_clicks"),
    prevent_initial_call=True
)
def return_manual(n_clicks):
    return dcc.send_file(
        "./data/goloco_user_manual.pdf"
        )
