import dash
from dash import html, callback
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash import dash_table
from dash.dash_table.Format import Format, Scheme
import pickle
import pandas as pd

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
    'width': '15%',
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
    'width': '15%',
    'borderWidth': '1px',
    'borderStyle': 'solid',
    'borderRadius': '1px',
}

df_dummy_2 = pd.DataFrame(columns=['ID', 'Source', 'Term ID', 'Term Name', 'p-value (adj)'])
df_dummy_3 = pd.DataFrame(columns=['HGNC ID', 'Ensembl ID'])

manhattan_overview = html.Div([
                html.P('A gene set enrichment analysis with gProfiler can be run by selecting genes from the z-score and volcano plot tabs in the "hits" page',
                style={'font-size': '1.2rem', 'font-style': 'italic'}),
                dbc.Row([
                dbc.Col([
                    html.Div([dbc.Card([dbc.FormGroup([
                    html.Br(),
                    dcc.Loading(id='manhattan-load', children=[dcc.Graph(id='manhattan-plot', figure = {"layout": {"height": 400, "width": 1000}})], type='default'),
                    html.Br(),
                    dcc.Loading(id = 'table_loading_2', children =[
                    html.Div(id='output-term-table', children = dash_table.DataTable(data=df_dummy_2.to_dict('records'),
                                        columns=[{'name': 'ID', 'id': 'ID', 'type':'text'},
                                                {'name': 'Source', 'id': 'Source', 'type': 'text'},
                                                {'name': 'Term ID', 'id': 'Term ID', 'type': 'text'},
                                                {'name': 'Term Name', 'id': 'Term Name', 'type': 'text'},
                                                {'name': 'p-value (adj)', 'id': 'p-value (adj)', 'type': 'numeric', 'format': Format(precision=2, scheme=Scheme.exponent)}],
                                        filter_action='native',
                                        sort_action='native',
                                        sort_mode='multi',
                                        style_table={
                                            'height': 300,
                                            'overflow': 'hidden',
                                            'textOverflow': 'ellipsis',
                                        },
                                        style_header={
                                                    'backgroundColor': 'lightgrey',
                                                    'font_family': 'Helvetica',
                                                    'font_size': '14px',
                                                    'paddingLeft': '20px',
                                                    'color': 'black',
                                                    'fontWeight': 'bold'
                                                },
                                        style_cell={'textAlign': 'left',                               
                                                    'font_family': 'Helvetica',
                                                    'font_size': '11px',
                                                    'paddingLeft': '20px',
                                                    },
                                        style_data={
                                            'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                                            'overflow': 'hidden',
                                            'tex4tOverflow': 'ellipsis',
                                            'minWidth': 50,
                                        }
                                        ), style={'padding': '20px'})], type = 'default'),
                ], style={'padding': '5px'})]),])
                ], width=9), 
                dbc.Col(dbc.Card([dbc.FormGroup([html.Div([
                    html.H4("Gene Mapping", className="display-5"),
                    html.P('Submit a gProfiler query to view gene mapping results:', style={'font-style': 'italic'}), 
                    html.Div(id='gprofile-gene-mapping',
                             children = [dcc.Loading(id='load_table', children=[dbc.Table.from_dataframe(df_dummy_3[['HGNC ID', 'Ensembl ID']], striped=True, bordered=True, hover=True, size='sm')], type='default'),], 
                                 style={
                                            "height": "725px",
                                            "overflow": "scroll",  # enable scrolling
                                        },),
                    html.Br(), 
                    ],
                         )])], style={'padding': '10px'}), width=3)]),
                         
                html.Br(),
                dbc.Row([
                    dbc.Col([
                        html.Div([
                            dbc.Button("+", id="collapse-button", outline=True, color="secondary", className="me-1", style={'margin-right': '10px'}, size="sm"),
                            dbc.Label('Expand to show query metadata', style={'font-size': '1.2rem'}),
                            html.Br(),
                        dbc.Collapse(
                            dbc.Card([dbc.FormGroup([html.Div(
                                children=[],
                                id="meta-data-info",
                                style={
                                    "width": "100%",
                                    "overflowX": "scroll",  # enable scrolling
                                }
                            )])])
                            ,
                            id="meta-collapse",
                            is_open=False,),
                        ])
                    ])
                ])

                ], style=SIDEBAR_STYLE)

layout = dbc.Container([
    html.Div([
    #dcc.Loading(id='warning_genes_load', children = [html.Div([html.Label(id='data_available_genes', style={'color':'#548235'})], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})], type='default'),
    html.Br(),
    dcc.Tabs(id="tabs-styled-with-inline-2", value='manhattan-overview', children = [
        dcc.Tab(label='''Overview''', value='manhattan-overview', style=tab_style, selected_style=tab_selected_style, children=manhattan_overview),
        #dcc.Tab(label='''Details''', value='linear-reg-tab', style=tab_style, selected_style=tab_selected_style, children=pairwise_gene_comp),
    ]),

    html.Div(id='dummy_div_enrichment'),],

    style={
        "transform": "scale(0.85)",
        "transform-origin": "top",
        "className": 'small'
    },
             ),

], fluid = True, style={'className': 'small'})

@callback(
    Output("meta-collapse", "is_open"),
    [Input("collapse-button", "n_clicks")],
    [State("meta-collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open