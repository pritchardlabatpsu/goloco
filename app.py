# imports
import os
import io
import base64
import json
import logging
import pandas as pd
import pickle
import numpy as np
import scipy.stats as stats
import plotly.express as px
import plotly.graph_objs as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
from urllib.parse import quote as urlquote
import boto3
import uuid
import base64
import dash
import dash_uploader as du
from dash import DiskcacheManager, CeleryManager, dash_table, ctx
from dash.dash_table.Format import Format
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from app_session import infer
import time
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.linear_model import LinearRegression
import random
from itertools import cycle
from collections import Counter
import tensorflow
import chronos
import seaborn as sns
import regex as re
import networkx as nx
from community import community_louvain
import phate
from dotenv import load_dotenv
from dash_extensions.enrich import callback, DashProxy, Output, Input, State, ServersideOutput, html, dcc, ServersideOutputTransform, MultiplexerTransform, FileSystemStore

load_dotenv()

# Directories
DepMap_DIRECTORY = './data'
L200_GENES = './data/l200_genes.pkl'
PREDICTION_DIRECTORY = './data/PC9SKMELCOLO_prediction.csv'
DepMap19Q4_DIRECTORY = './data/df_crispr_19q4.pkl'
DepMap19Q4_SAMPLE_INFO = './data/sample_info.pkl'
DepMap19Q4_GENES_DIRECTORY = './data/19q4_genes.csv'
DepMap19Q4_GENES_PICKLE = './data/19q4_genes.pkl'
DepMap19Q4_SUMMARY_DIRECTORY = './data/19q4_sum_stats.csv'
DepMap19Q4_ESSENTIAL_DIRECTORY = './data/common_essentials.csv'
DepMap19Q4_NONESSENTIAL_DIRECTORY = './data/nonessentials.csv'
DepMap19Q4_GENES_PICKLE = './data/19q4_genes.pkl'
AchillesCommonEssentials = './data/AchillesCommonEssentialControls.csv'
gene_stats_dir = './data/19q4_sum_stats.csv'

network_DIRECTORY = './data'
adjlist_DIRECTORY = './data/varExp_filtered.adjlist'
feat_sum_DIRECTORY = './data/feat_summary_varExp_filtered_class.csv'
modularity_DIRECTORY = './data/modularity.pkl'
cell_line_DIRECTORY = './data/cell_lines.pkl'

launch_uid = uuid.uuid4()

# AWS Bucket Information
s3_client = boto3.client('s3')
s3_bucket_name = 'ceres-infer'
s3 = boto3.resource('s3', aws_access_key_id=os.environ['AWS_ACCESS_KEY_ID'],
                          aws_secret_access_key=os.environ['AWS_SECRET_ACCESS_KEY'],
                          region_name='us-east-1')
my_bucket = s3.Bucket(s3_bucket_name)

# read depmap data from pickled directories
global dm_data
dm_data = pd.read_pickle(DepMap19Q4_DIRECTORY)
df_sample_data = pd.read_pickle(DepMap19Q4_SAMPLE_INFO)
l200_genes = pd.read_pickle(L200_GENES)

# read pickled data for biological network information and communities
global G
G = nx.read_adjlist(open(adjlist_DIRECTORY, 'rb'))

communities = community_louvain.best_partition(G)
nx.set_node_attributes(G, communities, 'modularity')

modularity = pickle.load(open(modularity_DIRECTORY, "rb"))

df_clusters = pd.read_csv('./data/feat_summary_varExp_filtered_class.csv')
landmark_genes = list(set(df_clusters['feature'].tolist()))   

# set-up redis and celery app; redis used as message broker for celery app backend
if 'REDIS_URL' in os.environ:
    logging.warning('Redis')
    # Use Redis & Celery if REDIS_URL set as an env variable
    from celery import Celery
    
    redis_url = os.environ['REDIS_URL']

    # set up celery app
    celery_app = Celery(__name__, 
                    BROKER_URL=redis_url, 
                    CELERY_RESULT_BACKEND=redis_url, 
                    BROKER_POOL_LIMIT=0, 
                    include=['app']
                    )

    # background callback manager used to queue tasks and return results
    background_callback_manager = CeleryManager(celery_app, cache_by=[lambda: launch_uid], expire=3600)
    my_backend = FileSystemStore(cache_dir='./cache')
    aws_available = True

else:
    # Diskcache for non-production apps when developing locally
    logging.warning('disk')

    import diskcache
    cache = diskcache.Cache()
    from celery import Celery

    redis_url = 'redis://localhost:6379/0'

    celery_app = Celery(__name__, 
                    BROKER_URL=redis_url, 
                    CELERY_RESULT_BACKEND=redis_url, 
                    BROKER_POOL_LIMIT=0, 
                    include=['app']
                    )

    background_callback_manager = CeleryManager(celery_app, cache_by=[lambda: launch_uid], expire=3600)
    my_backend = FileSystemStore(cache_dir='./cache')
    aws_available=False

# define dash app
app = DashProxy(__name__, use_pages=True, suppress_callback_exceptions=True,
                background_callback_manager=background_callback_manager,
                external_stylesheets=[dbc.themes.BOOTSTRAP],
                meta_tags=[{'name': 'viewport', 'content': 'width=2000px, initial-scale=0.8'}], 
                transforms=[ServersideOutputTransform(backend=my_backend), MultiplexerTransform()])

# define server for running app
server = app.server
#----------------------------------APP LAYOUT---------------------------------------------------------------------------------

app.layout = html.Div([
    dbc.NavbarSimple(
    children=[
        dbc.NavItem(dbc.NavLink("Home", href="/")),
        dbc.NavItem(dbc.NavLink("Predict", href="/predict")),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("Overview", href="/overview"),
                dbc.DropdownMenuItem("Genes", href="/genes"),
                dbc.DropdownMenuItem("Z-Scores", href="/zscore"),
            ],
            nav=True,
            in_navbar=True,
            label="Explore",
        ),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("Communities", href="/clusters"),
                dbc.DropdownMenuItem("Gene Walk", href="#"),
            ],
            nav=True,
            in_navbar=True,
            label="Pathways",
        ),
    ],
    brand="GO-LoCo",
    brand_href="/",
    color="primary",
    dark=True,
    style = {'padding-left': '50px',
             'padding-right': '50px'}
    ),

    dash.page_container,

    dcc.Store(id='L200-data-store', storage_type='session'),
    dcc.Store(id='load-df_pred-store', storage_type='session'),
    dcc.Store(id='input_filename', storage_type='session'),
    dcc.Store(id='total_hist-store', storage_type='memory'),
    dcc.Store(id='z_score_hist-store', storage_type='memory'),
    dcc.Loading(id='store-data-loading', children = [html.Div([dcc.Store(id='df_pred-store', storage_type='session')])], type='circle', fullscreen=True),
    dcc.Store(id='df_pred-store-tmp-1', storage_type='memory'),
    dcc.Store(id='df_pred-store-tmp-2', storage_type='memory'),
    dcc.Store(id='df_pred-store-tmp-3', storage_type='memory'),
    dcc.Store(id='experiment_labels', storage_type='session'),
    dcc.Store(id='df-counts-convert-tmp', storage_type='memory'),
    dcc.Store(id='df-counts-convert', storage_type='session'),
    dcc.Store(id='cell-line-or-experiment', storage_type='memory'),
    dcc.Store(id='cell-line-or-experiment-2', storage_type='memory'),
    #dcc.Store(id='gene_list', storage_type='session', data=gene_list),
    dcc.Location(id='url_name'),
    dcc.Download(id='download-pred'),
    dcc.Download(id='download-convert'),

    html.Div(id ='df_pred-store-html'),

    dbc.Modal(
        [
            dbc.ModalHeader("Congratulations"),
            dbc.ModalBody("Your genome wide inference has successfully completed. Select the button below to download your prediction. Navigate to the Explore tabs to visualize your results."),
            dbc.ModalFooter(
                dbc.Button(
                    "Download and Continue", id="completed", color="primary", className="ms-auto", n_clicks=0
                )
            ),
        ],
        id="completed-popup",
        is_open=False,
    ),

    dbc.Modal(
        [
            dbc.ModalHeader("Warning"),
            dbc.ModalBody("Are you sure you want to cancel your genome-wide inference?"),
            dbc.ModalFooter(
                dbc.Row([dbc.Col(dbc.Button(
                    "Nevermind", id="nevermind", color="primary", className="ms-auto", n_clicks=0)), 
                         dbc.Col(dbc.Button(
                    "Cancel", id="yes-cancel", color="danger", className="ms-auto", n_clicks=0
                ))]),
            ),
        ],
        id="cancel-popup",
        is_open=False,
    ),

    html.Div(id='dummy_div_app'), 

    dbc.Modal(
        [
            dbc.ModalHeader("Congratulations"),
            dbc.ModalBody("Converting your read counts to CHRONOS scores has completed. Please select the button below to download your data."),
            dbc.ModalFooter(
                dbc.Button(
                    "Download and Continue", id="conversion_completed", color="primary", className="ms-auto", n_clicks=0
                )
            ),
        ],
        id="conversion-completed-popup",
        is_open=False,
    ),
])

#-----------------------------------------GENERAL FUNCTIONS---------------------------------------------------------------

def save_file(name, content, directory):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(directory, name), "wb") as fp:
        fp.write(base64.decodebytes(data))
        
def uploaded_files(directory):
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(directory):
        path = os.path.join(directory, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files

def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "/download/{}".format(urlquote(filename))
    return html.A(filename, href=location)

def parse_contents(contents, filename):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return
    return df

def myround(x, base=5):
    return base * round(x/base)

def generate_uniprot_url(gene):
    return f'https://www.uniprot.org/uniprot/?query={gene}'

def create_gene_link(gene, website, name):
    if website == 'uniprot':
        url = generate_uniprot_url(gene)
    elif website == 'depmap':
        url = generate_depmap_url(gene)
    elif website == 'tumor_portal':
        url = generate_tumor_portal_url(gene)
    elif website == 'ncbi':
        url = generate_ncbi_url(gene)

    if name == 'gene':
        text = gene
    else:
        text = name
    
    return html.A(text, href=url, target='_blank')

def generate_depmap_url(gene):
    return f'https://depmap.org/portal/gene/{gene}?tab=overview'

def generate_tumor_portal_url(gene):
    return f'http://www.tumorportal.org/view?geneSymbol={gene}'

def generate_ncbi_url(gene):
    return f'https://www.ncbi.nlm.nih.gov/gene/?term={gene}'

def render_link(value):
    return dcc.Link(value, href='https://www.uniprot.org/uniprot/'+value)

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    c = colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])
    c = mc.to_hex(c)
    return c

def adjust_lightness_2(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    c = colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2], 0.5)
    return c

def find_alpha_05(data):
    # Calculate the sample mean and standard deviation
    sample_mean = np.mean(data)
    sample_std = np.std(data, ddof=1)

    # Calculate the critical values for a two-tailed test with alpha = 0.05
    alpha = 0.05
    z_crit = stats.norm.ppf(alpha / 2), stats.norm.ppf(1 - alpha / 2)

    # Calculate the corresponding confidence interval
    conf_interval = sample_mean + np.array(z_crit) * sample_std
    return conf_interval

def check_dna(list_a):
    """
    Check if all columns in a pandas dataframe only contain strings and in those strings only contain characters 'T', 'C', 'G', 'A' or 'U'
    
    Parameters:
    df (pd.DataFrame): DataFrame to be checked
    
    Returns:
    bool: True if all columns contain only valid characters, False otherwise.
    """
    valid_chars = set('TCGAU')
    for col in list_a:
        if not isinstance(col, str) and set(col.upper()).issubset(valid_chars):
            return False
    return True

def check_int(df):
    """
    Check if all values in a pandas dataframe are integers
    
    Parameters:
    df (pd.DataFrame): DataFrame to be checked
    
    Returns:
    bool: True if all values are integers, False otherwise.
    """
    return df.applymap(lambda x: isinstance(x, int)).all().all()

def check_format(df, column_name):
    """
    Checks if every value in a column of a pandas dataframe follows the format "geneSym (geneID)"
    
    Parameters:
    df (pandas.DataFrame): The pandas dataframe
    column_name (str): The name of the column to check
    
    Returns:
    bool: True if every value in the column follows the format, False otherwise
    """
    pattern = r'^[a-zA-Z0-9]+ \(\d+\)$'  # regular expression pattern for the format
    
    # Check each value in the column against the pattern
    for value in df[column_name]:
        if not re.match(pattern, value):
            return False
    
    return True

def calc_z_score(stats_data, gene2analyz, pred_score):

    x_avg = stats_data.loc[gene2analyz]['avg']

    x_std = stats_data.loc[gene2analyz]['std']

    z_score = (pred_score - x_avg) / x_std
    return z_score

def return_s3_object(s3, s3_bucket_name, pkl_filename):
    obj = s3.Object(s3_bucket_name, pkl_filename)
    data = obj.get()['Body'].read()
    pkl_data = io.BytesIO(data)
    return pkl_data
#----------------------------------------CALLBACKS----------------------------------------------------------------

@app.callback(
    Output('warning-popup', 'is_open'),
    [Input('submit-inference', 'n_clicks'), Input('close', 'n_clicks')],
    [State('upload-L200-data', 'filename'), State('warning-popup', 'is_open')],
    prevent_initial_call = True
)
def toggle_modal(n1, n2, L200_filename, is_open):
    if L200_filename is None:
        if n1 or n2:
            return not is_open
        return is_open

@app.callback(
    Output('continue-popup', 'is_open'),
    [Input('submit-inference', 'n_clicks'), Input('continue', 'n_clicks')],
    [State('upload-L200-data', 'filename'), State('continue-popup', 'is_open')],
    prevent_initial_call = True
)
def toggle_modal(n1, n2, L200_filename, is_open):
    if L200_filename is not None:
        if n1 or n2:
            return not is_open
        return is_open

@app.callback(
    Output('cancel-popup', 'is_open'),
    [Input('nevermind', 'n_clicks'), Input('yes-cancel', 'n_clicks'), Input('cancel-inference', 'n_clicks')],
    [State('cancel-popup', 'is_open')],
    prevent_initial_call = True)
def toggle_modal(n1, n2, n3, is_open):
    if n1 or n2 or n3:
        return not is_open
    return is_open

@app.callback(
    Output('warning-popup-2', 'is_open'),
    [Input('submit-prediction', 'n_clicks'), Input('close_2', 'n_clicks')],
    [State('upload-df_pred-data', 'filename'), State('warning-popup-2', 'is_open')],
    prevent_initial_call = True
)
def toggle_modal(n1, n2, L200_filename, is_open):
    if L200_filename is None:
        if n1 or n2:
            return not is_open
        return is_open

@app.callback(Output('completed-popup', 'is_open'),
              Input('df_pred-store', 'data'),
              Input('completed', 'n_clicks'),
              State('completed-popup', 'is_open'),
              prevent_initial_call = True)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open

@app.callback(Output('conversion-completed-popup', 'is_open'),
              Input('df-counts-convert', 'data'),
              Input('conversion_completed', 'n_clicks'),
              State('conversion-completed-popup', 'is_open'),
              prevent_initial_call = True)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open

#---------------------------------------------------------------PREDICT PAGE CALLBACKS-------------------------------------------------------

@app.callback(
    Output('L200-data-store', 'data'),
    Output('file_upload_label', 'children'),
    Output('submit-inference', 'disabled'),
    Input('upload-L200-data', 'contents'),
    State('upload-L200-data', 'filename'),
    prevent_initial_call = True)
def store_l200_data(l200_contents, l200_filename):
    if l200_contents is None or l200_filename is None:
        raise PreventUpdate()
    else:
        try:
            data = parse_contents(l200_contents, l200_filename)
            df_l200 = pd.DataFrame(data.to_dict('records'))
            genes = df_l200['feature'].values.tolist()
            df_l200 = df_l200.set_index('feature')
            exp = df_l200.columns.tolist()
            na = df_l200.isnull().values.any()
            if len(genes)==200 and Counter(genes)==Counter(l200_genes) and len(exp)>0 and not na:
                return data.to_dict('records'), [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-check-circle-fill me-2"),
                                "File Uploaded Succesfully: " + l200_filename,
                            ],
                            color="success",
                            className="d-flex align-items-center",
                        ),], False
            else:
                return None, [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-x-octagon-fill me-2"),
                                "Processing Error: please review formatting of " + l200_filename + " or try another file",
                            ],
                            color="danger",
                            className="d-flex align-items-center",
                        ),], True
        except Exception as err:
            return None, [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-x-octagon-fill me-2"),
                                "Processing Error: please review formatting of " + l200_filename + " or try another file",
                            ],
                            color="danger",
                            className="d-flex align-items-center",
                        ),], True

@app.callback(
    Output('pred_upload_label', 'children'),
    Output('submit-prediction', 'disabled'),
    Input('upload-df_pred-data', 'contents'),
    State('upload-df_pred-data', 'filename'),
    prevent_initial_call = True)
def store_l200_data(df_pred_contents, df_pred_filename):
    if df_pred_contents is None or df_pred_filename is None:
        raise PreventUpdate()
    else:
        try:
            data = parse_contents(df_pred_contents, df_pred_filename)
            x = data.to_dict('dict')
            experiments = list(x.keys())
            exp = [i.replace(' (CERES Pred)', '') for i in experiments if '(CERES Pred)' in i]
            zscore = [i.replace(' (Z-Score)', '') for i in experiments if '(Z-Score)' in i]
            cols = [i for i in experiments if i in ['gene_category', 'gene', 'avg', 'std']]
            if len(exp) > 0 and len(zscore) > 0 and len(cols)>0:
                return [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-check-circle-fill me-2"),
                                "File Uploaded Succesfully: " + df_pred_filename,
                            ],
                            color="success",
                            className="d-flex align-items-center",
                        ),], False
            else: 
                return [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-x-octagon-fill me-2"),
                                "Processing Error: please review formatting of " + df_pred_filename + " or try another file",
                            ],
                            color="danger",
                            className="d-flex align-items-center",
                        ),], True
        except Exception as err:
            return [html.Br(),
            dbc.Alert(
                [
                    html.I(className="bi bi-x-octagon-fill me-2"),
                    "Processing Error: please review formatting of " + df_pred_filename + " or try another file",
                ],
                color="danger",
                className="d-flex align-items-center",
            ),], True

@app.callback(
    Output('upload_reads_label', 'children'),
    Input('upload-reads-data', 'contents'),
    State('upload-reads-data', 'filename'),
    prevent_initial_call = True)
def readcounts_data_upload(read_data, read_filename):
    if read_data is None or read_filename is None:
        raise PreventUpdate()
    else:
        try:
            warning_messages = []
            readcounts = parse_contents(read_data, read_filename)
            k = readcounts.columns[0] == 'Unnamed: 0'
            if not k:
                warning_messages.append('First column does not appear to be an unnamed index')
            readcounts = readcounts.set_index(readcounts.iloc[:, 0])
            del readcounts[readcounts.columns[0]]
            a = check_dna(readcounts.columns.tolist())
            if not a:
                warning_messages.append('Column names do not match expected formats of sgRNA sequences')
            b = check_int(readcounts)
            if not b:
                warning_messages.append('At least one readcount value missing or not an integer')

            if a and b and k:
                return [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-check-circle-fill me-2"),
                                "Everything looks good! File Uploaded Succesfully: " + read_filename,
                            ],
                            color="success",
                            className="d-flex align-items-center",
                        ),]
            else:
                return [html.Br(),
                        dbc.Alert(
                            [
                                html.Div([
                                html.I(className="bi bi-x-octagon-fill me-2"),
                                "Processing Error: please review formatting of " + read_filename + " or try another file",
                                html.Br(),
                                dbc.ListGroup(
                                [
                                    dbc.ListGroupItem([message]) for message in warning_messages
                                        ],),])
                                ],
                            color="danger",
                            className="d-flex align-items-center",
                        ),]
        except Exception as err:
            return [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-x-octagon-fill me-2"),
                                "Processing Error: please review formatting of " + read_filename + " or try another file",
                            ],
                            color="danger",
                            className="d-flex align-items-center",
                        ),]

@app.callback(
    Output('upload_sequence_map_label', 'children'),
    Input('upload-sequence-map-data', 'contents'),
    State('upload-sequence-map-data', 'filename'),
    prevent_initial_call = True)
def readcounts_data_upload(sequence_data, sequence_filename):
    if sequence_data is None or sequence_filename is None:
        raise PreventUpdate()
    else:
        try:
            warning_messages = []
            sequencemap = parse_contents(sequence_data, sequence_filename)
            column_list = ['sequence_ID', 'cell_line_name', 'pDNA_batch', 'days']
            a = Counter(sequencemap.columns.tolist()) == Counter(column_list)
            if not a:
                warning_messages.append('At least one required column with correct header is missing')
            sequencemap_small = sequencemap[sequencemap['cell_line_name']=='pDNA']
            batches = set(sequencemap.pDNA_batch.values.tolist())
            n_batches = set(sequencemap_small.pDNA_batch.values.tolist())
            b = Counter(batches)==Counter(n_batches)
            if not b:
                warning_messages.append('At least one pDNA batch is missing a pDNA measurment')
            
            if a and b:
                return [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-check-circle-fill me-2"),
                                "Everything looks good! File Uploaded Succesfully: " + sequence_filename,
                            ],
                            color="success",
                            className="d-flex align-items-center",
                        ),]
            else:
                return [html.Br(),
                        dbc.Alert(
                            [
                                html.Div([
                                html.I(className="bi bi-x-octagon-fill me-2"),
                                "Processing Error: please review formatting of " + sequence_filename + " or try another file",
                                html.Br(),
                                dbc.ListGroup(
                                [
                                    dbc.ListGroupItem([message]) for message in warning_messages
                                        ],)])
                            ],
                            color="danger",
                            className="d-flex align-items-center",
                        ),]
        except Exception as err:
            return [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-x-octagon-fill me-2"),
                                "Processing Error: please review formatting of " + sequence_filename + " or try another file",
                            ],
                            color="danger",
                            className="d-flex align-items-center",
                        ),]

@app.callback(
    Output('upload_guide_map_label', 'children'),
    Input('upload-guide-map-data', 'contents'),
    State('upload-guide-map-data', 'filename'),
    prevent_initial_call = True)
def readcounts_data_upload(guide_data, guide_filename):
    if guide_data is None or guide_filename is None:
        raise PreventUpdate()
    else:
        try:
            warning_messages = []
            guidemap = parse_contents(guide_data, guide_filename)
            column_list = ['sgrna', 'gene']
            a = set(guidemap.columns.tolist()).issubset(column_list)
            if not a:
                warning_messages.append('At least one required column with correct header is missing')
            b = len(set(guidemap.sgrna.values.tolist())) == len(guidemap.sgrna.values.tolist())
            if not b:
                warning_messages.append('At least one sgRNA sequence is repeated or mapped twice to different genes')
            c = check_format(guidemap, 'gene')
            if not c:
                warning_messages.append('At least one gene does not follow the correct format "geneSym (geneID)"')
            d = check_dna(guidemap['sgrna'].values.tolist())
            if not d:
                warning_messages.append('At least one sgrna sequence does not match expected format of sgRNA sequences')

            if a and b and c and d:
                return [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-check-circle-fill me-2"),
                                "Everything looks good! File Uploaded Succesfully: " + guide_filename,
                            ],
                            color="success",
                            className="d-flex align-items-center",
                        ),]
            else:
                return [html.Br(),
                        dbc.Alert(
                            [
                                html.Div([
                                html.I(className="bi bi-x-octagon-fill me-2"),
                                "Processing Error: please review formatting of " + guide_filename + " or try another file",
                                html.Br(),
                                dbc.ListGroup(
                                [
                                    dbc.ListGroupItem([message]) for message in warning_messages
                                        ],)])
                            ],
                            color="danger",
                            className="d-flex align-items-center",
                        ),]
        except Exception as err:
            return [html.Br(),
                        dbc.Alert(
                            [
                                html.I(className="bi bi-x-octagon-fill me-2"),
                                "Processing Error: please review formatting of " + guide_filename + " or try another file",
                            ],
                            color="danger",
                            className="d-flex align-items-center",
                        ),]

@app.callback(
    Output('upload_all_convert_files', 'children'),
    Output('submit-convert', 'disabled'),
    Input('upload-reads-data', 'contents'),
    Input('upload-reads-data', 'filename'),
    Input('upload-sequence-map-data', 'contents'),
    Input('upload-sequence-map-data', 'filename'),
    Input('upload-guide-map-data', 'contents'),
    Input('upload-guide-map-data', 'filename'),)
def upload_all(read_data, read_filename, sequence_data, sequence_filename, guide_data, guide_filename):
    if read_data is None or sequence_data is None or guide_data is None:
        return dash.no_update
    else:
        readcounts = parse_contents(read_data, read_filename)
        k = readcounts.columns[0] == 'Unnamed: 0'
        readcounts = readcounts.set_index(readcounts.iloc[:, 0])
        del readcounts[readcounts.columns[0]]
        a = check_dna(readcounts.columns.tolist())
        b = check_int(readcounts)

        sequencemap = parse_contents(sequence_data, sequence_filename)
        column_list = ['sequence_ID', 'cell_line_name', 'pDNA_batch', 'days']
        c = Counter(sequencemap.columns.tolist()) == Counter(column_list)
        sequencemap_small = sequencemap[sequencemap['cell_line_name']=='pDNA']
        batches = set(sequencemap.pDNA_batch.values.tolist())
        n_batches = set(sequencemap_small.pDNA_batch.values.tolist())
        d = Counter(batches)==Counter(n_batches)

        guidemap = parse_contents(guide_data, guide_filename)
        column_list = ['sgrna', 'gene']
        e = set(guidemap.columns.tolist()).issubset(column_list)
        f = len(set(guidemap.sgrna.values.tolist())) == len(guidemap.sgrna.values.tolist())
        g = check_format(guidemap, 'gene')
        h = check_dna(guidemap['sgrna'].values.tolist())

        warning_messages = []

        i = Counter(readcounts.columns.tolist())==Counter(guidemap['sgrna'].values.tolist())
        j = Counter(readcounts.index.values.tolist())==Counter(sequencemap['sequence_ID'].values.tolist())

        if not i:
            warning_messages.append('Columns in read counts dataframe do not match sgRNA reads in guide map dataframe')
        if not j:
            warning_messages.append('Indexes in read counts dataframe do not match sequence_IDs in sequence map dataframe')

        if a and b and c and d and e and f and g and h and k:
            if i and j:
                return [html.Br(),
                            dbc.Alert(
                                [
                                    html.I(className="bi bi-check-circle-fill me-2"),
                                    "Everything looks good! Ready to convert scores!",
                                ],
                                color="success",
                                className="d-flex align-items-center",
                            ),], False
            else:
                return [dbc.Alert(
                            [
                            html.Div([
                            dbc.ListGroup(
                            [
                                dbc.ListGroupItem([message]) for message in warning_messages
                                    ],),
                             ],),
                                                ],
                             color="danger",
                             className="d-flex align-items-center",
                             style={'align-items': 'center'}
                       ),], True
        else:
            if i and j:
                return [dbc.Alert(
                            [
                            html.Div([
                            html.I(className="bi bi-info-circle-fill me-2"),
                            'Something does not look right, check your dataframes'
                             ],),
                                                ],
                             color="danger",
                             className="d-flex align-items-center",
                             style={'align-items': 'center'}
                       ),], True
            else:
                return [dbc.Alert(
                            [
                            html.Div([
                            html.I(className="bi bi-info-circle-fill me-2"),
                            html.Br(),
                            dbc.ListGroup(
                            [
                                dbc.ListGroupItem([message]) for message in warning_messages
                                    ],),
                             ],),
                                                ],
                             color="danger",
                             className="d-flex align-items-center",
                             style={'align-items': 'center'}
                       ),], True

@app.callback(
    Output('df-counts-convert-tmp', 'data'),
    Output('progress_bar_convert', 'children'),
    Input('submit-convert', 'n_clicks'),
    State('upload-reads-data', 'contents'),
    State('upload-reads-data', 'filename'),
    State('upload-sequence-map-data', 'contents'),
    State('upload-sequence-map-data', 'filename'),
    State('upload-guide-map-data', 'contents'),
    State('upload-guide-map-data', 'filename'),
    running = [(Output('submit-convert', 'disabled'), True, False)],
    background = True,
    prevent_initial_call = True
    )
def run_convert(n_clicks, read_data, read_filename, sequence_data, sequence_filename, guide_data, guide_filename):
    if not n_clicks:
        return dash.no_update
    elif read_data is None or sequence_data is None or guide_data is None:
        return dash.no_update
    elif read_data is not None and sequence_data is not None and guide_data is not None:
        readcounts = parse_contents(read_data, read_filename)
        guidemap = parse_contents(guide_data, guide_filename)
        sequencemap = parse_contents(sequence_data, sequence_filename)

        readcounts = readcounts.set_index(readcounts.iloc[:, 0])
        del readcounts[readcounts.columns[0]]

        chronos.nan_outgrowths(readcounts=readcounts, guide_gene_map=guidemap, sequence_map=sequencemap)
        model = chronos.Chronos(
                sequence_map={'brunello': sequencemap},
                guide_gene_map={'brunello': guidemap},
                readcounts={'brunello': readcounts}
            )

        model.train(801)
        gene_effects = model.gene_effect

        if len(set(guidemap.gene.values.tolist())) == 200:
            cell_lines = gene_effects.index.tolist()
            gene_effects = gene_effects.T.reset_index()
            gene_effects = gene_effects.rename(columns={'index': 'gene'})
            gene_effects['gene_ID'] = gene_effects['gene'].str.extract(r'\((.*)\)')[0]
            gene_effects['gene_sym'] = gene_effects['gene'].str.split(pat=' ', n=1, expand=True)[0]

            key = 'Gene'

            essntl_22q4 = pd.read_csv(AchillesCommonEssentials)
            essntl_22q4[key] = essntl_22q4[key].str.extract(r'\((.*)\)')[0]
            essntl_22q4['gene_type'] = 'essential'

            gene_effects = gene_effects.merge(essntl_22q4, how='left', right_on=key, left_on='gene_ID')
            gene_effects = gene_effects.drop([key], axis=1)
            gene_effects['gene_type'] = gene_effects['gene_type'].fillna('conditional')

            ess_mean = -0.974066437522624
            cond_ess_mean = -0.2591980795361984

            for line in cell_lines:
                med_cond_ess_L200 = gene_effects.loc[gene_effects['gene_type']=='conditional', line].median()
                med_essntl_L200 = gene_effects.loc[gene_effects['gene_type']=='essential', line].median()
    
                m = (cond_ess_mean - ess_mean) / (med_cond_ess_L200 - med_essntl_L200)
                b = cond_ess_mean - m * med_cond_ess_L200
    
                gene_effects[line] = m * gene_effects[line] + b

            gene_effects = gene_effects[['gene'] + cell_lines]
            gene_effects = gene_effects.set_index('gene')

        user_id = str(uuid.uuid1())
        pkl_filename = "chronos/" + user_id + "_l200_gene_effect_conversion.pkl"
        pickle_buffer = io.BytesIO()
        gene_effects.to_pickle(pickle_buffer)
        s3.Object(s3_bucket_name, pkl_filename).put(Body=pickle_buffer.getvalue())

        return pkl_filename, None

@app.callback(ServersideOutput('df-counts-convert', 'data'),
              Input('df-counts-convert-tmp', 'data'),)
def update_df_convert(pkl_filename):
    if pkl_filename is None:
        return dash.no_update
    else:
        obj = s3.Object(s3_bucket_name, pkl_filename)
        data = obj.get()['Body'].read()
        gene_effects = pickle.load(io.BytesIO(data))                 
        return gene_effects.to_dict('dict')

@app.callback(Output('download-convert', 'data'),
              Input('conversion_completed', 'n_clicks'),
              State('df-counts-convert', 'data'),
              prevent_initial_call = True)
def update_table(n_clicks, convert_data):
    gene_effects = pd.DataFrame(convert_data)
    return dcc.send_data_frame(gene_effects.to_csv, 'CHRONOS_Scores.csv')

@app.callback(
    Output('df_pred-store-tmp-3', 'data'),
    Output('experiment_labels', 'data'),
    Output('input_filename', 'data'),
    Input('continue', 'n_clicks'),
    State('upload-L200-data', 'filename'),
    State('L200-data-store', 'data'),
    running = [(Output('submit-inference', 'disabled'), True, False),
               (Output('load-prediction', 'disabled'), True, False),
               (Output('submit-prediction', 'disabled'), True, False),
               (Output('eta_label', 'hidden'), False, True),
               (Output("progress_bar", "style"), {"visibility": "visible", "height": "20px"}, {"visibility": "hidden", "height": "20px"},),
               (Output("cancel-inference", "style"), {"visibility": "visible"}, {"visibility": "hidden"},),
               (Output("cancel-inference", 'disabled'), False, True)
               ],
    cancel = [Input("yes-cancel", "n_clicks")],
    progress = [Output("progress_bar", "value"), Output("progress_bar", "label"), Output("eta_label", "children")],
    background = True, # for deployment
    prevent_initial_call = True)
def run_inference(set_progress, con_n_clicks, L200_filename, L200_data): #set_progress, 
    if not con_n_clicks:
        return dash.no_update
    elif L200_filename is None or L200_data is None:
        return dash.no_update
    elif L200_filename is not None and L200_data is not None:
        set_progress((0, "0 %", "Estimated Time Remaining: " + "inf"))
        scope = pd.read_pickle('./data/19q4_genes.pkl')
        df_l200 = pd.DataFrame(L200_data)
        a = infer(df_l200)
        genes = []
        ceres_pred = np.zeros(shape=(len(scope),len(a.experiments)))
        logging.warning(ceres_pred.shape)
        z_scores = np.zeros(shape=(len(scope),len(a.experiments)))
        logging.warning(z_scores.shape)
        x_avgs = []
        x_stds = []
        t1 = time.perf_counter()
        for i, gene in enumerate(scope):
            logging.warning(gene)
            pred = a.infer_gene(gene, aws_s3=True)
            x_avg, x_std, z_score = a.calc_z_score(gene, pred)
            genes.append(gene)
            ceres_pred[i] = pred
            z_scores[i] = z_score
            x_avgs.append(x_avg)
            x_stds.append(x_std)
            p = 100 * (i+1) / len(scope)
            p = np.round(p) - 1
            t2 = time.perf_counter()
            t_avg = (t2 - t1) / (i + 1)
            seconds = t_avg * (len(scope) - i)
            eta = time.strftime("%H:%M:%S", time.gmtime(seconds))
            set_progress((p, f'{p} %', "Estimated Time Remaining: " + eta))
        df_pred = pd.DataFrame()
        set_progress((99, "100 %", "Estimated Time Remaining: " + "00:00:00"))
        df_pred['gene'] = genes
        df_pred = pd.merge(df_pred, a.gene_categories, on='gene', how='left')
        df_pred['gene_category'] = df_pred['gene_category'].replace(np.nan, 'conditional essential')
        df_pred['avg'] = x_avgs
        df_pred['std'] = x_stds
        for i, exp in enumerate(a.experiments):
            df_pred[exp + ' (CERES Pred)'] = ceres_pred[:, i]
            df_pred[exp + ' (Z-Score)'] = z_scores[:, i]
        df_pred.set_index('gene')
        user_id = str(uuid.uuid1())
        pkl_filename = "predictions/" + user_id + "_genome_wide_prediction.pkl"
        #df_pred = pd.read_csv(os.path.join(OUTPUT_DIRECTORY, 'prediction.csv')) # for testing
        pickle_buffer = io.BytesIO()
        df_pred.to_pickle(pickle_buffer)
        s3.Object(s3_bucket_name, pkl_filename).put(Body=pickle_buffer.getvalue())
        logging.warning(pkl_filename)
        set_progress((100, "100 %", "Estimated Time Remaining: " + "00:00:00"))
        return pkl_filename, a.experiments, L200_filename

@app.callback(ServersideOutput('df_pred-store-tmp-1', 'data'),
              Output('experiment_labels', 'data'),
              Output('input_filename', 'data'),
              Input('load-prediction', 'n_clicks'),
              prevent_initial_call=True)
def load_prediction(n_clicks):
    if not n_clicks:
        logging.warning('rasing prevent update')
        raise PreventUpdate()
    else:
        df_pred = pd.read_csv(PREDICTION_DIRECTORY)
        experiments = df_pred.columns.tolist()
        experiments = [i.replace(' (CERES Pred)', '') for i in experiments if '(CERES Pred)' in i]
        logging.warning(experiments)
        return df_pred.to_dict('dict'), experiments, 'preloaded'

@app.callback(ServersideOutput('df_pred-store-tmp-2', 'data'),
              Output('experiment_labels', 'data'),
              Output('input_filename', 'data'),
              Input('submit-prediction', 'n_clicks'),
              State('upload-df_pred-data', 'contents'),
              State('upload-df_pred-data', 'filename'),
              prevent_initial_call = True)
def submit_prediction(n_clicks, df_pred_contents, df_pred_filename):
    if not n_clicks or df_pred_contents is None or df_pred_filename is None:
        logging.warning('raising prevent update')
        raise PreventUpdate()
    else:
        data = parse_contents(df_pred_contents, df_pred_filename)
        x = data.to_dict('dict')
        experiments = list(x.keys())
        experiments = [i.replace(' (CERES Pred)', '') for i in experiments if '(CERES Pred)' in i]
        logging.warning(experiments)
        return x, experiments, df_pred_filename

@app.callback(ServersideOutput('df_pred-store', 'data'),
              Input('df_pred-store-tmp-1', 'data'),
              Input('df_pred-store-tmp-2', 'data'),
              Input('df_pred-store-tmp-3', 'data'))
def update_df_pred_store(data1, data2, pkl_filename):
    trigger_id = ctx.triggered_id
    if trigger_id=='df_pred-store-tmp-1':
        return data1
    elif trigger_id=='df_pred-store-tmp-2':
        return data2
    elif trigger_id=='df_pred-store-tmp-3':
        logging.warning('****')
        obj = s3.Object(s3_bucket_name, pkl_filename)
        data = obj.get()['Body'].read()
        df_pred = pickle.load(io.BytesIO(data))                 
        return df_pred.to_dict('dict')

@app.callback(Output('download-pred', 'data'),
              Input('completed', 'n_clicks'),
              State('df_pred-store', 'data'),
              State('input_filename', 'data'),
              prevent_initial_call = True)
def update_table(n_clicks, pred_data, filename):
    df_pred = pd.DataFrame(pred_data)
    filename = filename.replace('.csv', '')
    filename = filename.replace('.xlsx', '')
    filename = filename.replace('.xls', '')
    filename = filename + '_prediction.csv'
    return dcc.send_data_frame(df_pred.to_csv, filename)

#---------------------------------------------------------- OVERVIEW PAGE CALLBACKS---------------------------------------------------------

@app.callback(Output('output-inference-table', 'children'),
              Input('choose_experiment_1', 'value'),
              State('df_pred-store', 'data'),
              prevent_initial_call = True)
def update_table(experiment, pred_data):
    if experiment is None or pred_data is None:
        return dash.no_update
    df_pred = pd.DataFrame(pred_data)
    columns = ['gene', 'gene_category', 'avg', 'std'] + [e + ' (CERES Pred)' for e in experiment] + [e + ' (Z-Score)' for e in experiment]
    df_pred = df_pred[columns]
    all_colors = {}
    all_colors.update(mcolors.TABLEAU_COLORS)
    ceres_ids = [e + ' (CERES Pred)' for e in experiment]
    z_score_ids = [e + ' (Z-Score)' for e in experiment]
    color_dict = dict(zip(["conditional essential", "common nonessential", "common essential"], mcolors.TABLEAU_COLORS))

    column_list = [{'name': 'Gene', 'id': 'gene', 'presentation': 'markdown'},
                   {'name': 'Gene Category', 'id': 'gene_category', 'type': 'text'},
                   {'name': 'Average', 'id': 'avg', 'type': 'numeric', 'format':Format(precision=2),},
                   {'name': 'Standard Deviation', 'id': 'std', 'type': 'numeric', 'format':Format(precision=2),}]

    for e in experiment:
        a_dict = {}
        a_dict['name'] = e + '\n (CERES)'
        a_dict['id'] = e + ' (CERES Pred)'
        a_dict['type'] = 'numeric'
        a_dict['format'] = Format(precision=2)

        column_list.append(a_dict)

        a_dict = {}
        a_dict['name'] = e + '\n (Z-Score)'
        a_dict['id'] = e + ' (Z-Score)'
        a_dict['type'] = 'numeric'
        a_dict['format'] = Format(precision=2)

        column_list.append(a_dict)

    table = dash_table.DataTable(data=df_pred.to_dict('records'),
                                              columns=column_list,
                                              filter_action='native',
                                              sort_action='native',
                                              sort_mode='multi',
                                              style_cell={'textAlign': 'left',                               
                                                          'font_family': 'Arial',
                                                          'font_size': '16px',
                                                          'paddingLeft': '20px',
                                                          'padding': '5px',
                                                          'padding-left': '20px',},
                                              style_header={
                                                            'backgroundColor': 'dodgerblue',
                                                            'font_family': 'sans-serif',
                                                            'font_size': '18px',
                                                            'paddingLeft': '20px',
                                                            'padding': '5px',
                                                            'color': 'white',
                                                            'fontWeight': 'bold',
                                                            'padding-left': '20px',
                                                            'height': 'auto',
                                                            'whiteSpace': 'normal',
                                                            'maxWidth': '100px',
                                                        },
                                              style_table={
                                                    'height': 800,
                                                    'overflow': 'hidden',
                                                    'textOverflow': 'ellipsis',
                                              },
                                              style_data={
                                                    'width': '50px', 'minWidth': '50px', 'maxWidth': '150px',
                                                    'overflow': 'hidden',
                                                    'tex4tOverflow': 'ellipsis',
                                                    'minWidth': 50,
                                                    'whiteSpace': 'normal',
                                                    'height': 'auto',
                                                    'lineHeight': '15px'
                                              },
                                              style_data_conditional=([{'if': {
                                                                              'filter_query': '{gene_category} = "common essential"',
                                                                              'column_id': 'gene_category'
                                                                               },
                                                                               'backgroundColor': adjust_lightness(all_colors[color_dict["common essential"]], 1.7),
                                                                               'color': 'white',
                                                                               'fontWeight': 'bold'},
                                                                       {'if': {
                                                                              'filter_query': '{gene_category} = "common nonessential"',
                                                                              'column_id': 'gene_category'
                                                                               },
                                                                               'backgroundColor': adjust_lightness(all_colors[color_dict["common nonessential"]], 1.3),
                                                                               'color': 'white',
                                                                               'fontWeight': 'bold'},
                                                                       {'if': {
                                                                              'filter_query': '{gene_category} = "conditional essential"',
                                                                              'column_id': 'gene_category'
                                                                               },
                                                                               'backgroundColor': adjust_lightness(all_colors[color_dict["conditional essential"]], 1.7),
                                                                               'color': 'white',
                                                                               'fontWeight': 'bold'},
                                                                           ]),
                                              style_as_list_view=True,
                                              markdown_options={"html": True},
                                              )

    for ceres_id in ceres_ids:
        rule = {
            'if': {
                'filter_query': f'{{{ceres_id}}} < -1',
                'column_id': ceres_id,
            },
            'backgroundColor': adjust_lightness('#FF4136', 1.4),
            'color': 'white',
            'fontWeight': 'bold',
        }
        table.style_data_conditional.append(rule)

    for z_score_id in z_score_ids:
        rule = {
            'if': {
                'filter_query': f'{{{z_score_id}}} < -1',
                'column_id': z_score_id,
            },
            'backgroundColor': adjust_lightness('#FF4136', 1.4),
            'color': 'white',
            'fontWeight': 'bold',
        }
        table.style_data_conditional.append(rule)

    return html.Div([table])

@app.callback(Output('total_dist_graph','figure'),
              Output('z_score_dist_graph','figure'),
              Output('total_dist_graph_pdf','figure'),
              Output('z_score_dist_graph_pdf','figure'),
              Input('choose_experiment_2', 'value'),
              Input('gene_cat_checklist_1', 'value'),
              State('df_pred-store', 'data'),
              prevent_initial_call = True)
def update_overview(experiment, categories, pred_data):
    if experiment is None or categories is None or pred_data is None:
        return dash.no_update
    df_pred = pd.DataFrame(pred_data)
    all_colors = {}
    all_colors.update(mcolors.TABLEAU_COLORS)
    color_dict = dict(zip(["conditional essential", "common nonessential", "common essential"], mcolors.TABLEAU_COLORS))
    
    fig = go.Figure()

    fig.update_layout(
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'))
    fig.update_layout(plot_bgcolor='white')
    fig.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5)
    ymax = 0

    for i, cat in enumerate(categories):
        hist_graph = go.Histogram(x=df_pred[df_pred['gene_category']==cat][experiment + ' (CERES Pred)'].values, 
                                  name=cat, 
                                  marker_color=all_colors[color_dict[cat]],
                                  opacity=0.50)
        fig.add_trace(hist_graph)
        bins = fig.full_figure_for_development().data[0].xbins
        nbins = round((bins['end'] - bins['start'])/bins['size'])
        counts, edges = np.histogram(df_pred[df_pred['gene_category']==cat][experiment + ' (CERES Pred)'].values, bins=nbins)
        y = max(counts)
        if y > ymax:
            ymax = y
        ymax += 10
    for cat in categories:
        x_avg = np.average(df_pred[df_pred['gene_category']==cat][experiment + ' (CERES Pred)'].values)
        fig.add_trace(go.Scatter(x=[x_avg, x_avg], 
                                 y = [0, ymax], 
                                 mode='lines',
                                 line=dict(color=adjust_lightness(all_colors[color_dict[cat]], 0.9), width=2, dash='dot'),
                                 name=cat + ' (average)'))
    fig.update_layout(barmode='overlay')
    fig.update_layout(
    xaxis_title="Predicted CERES Score",
    yaxis_title="Number of Genes",
    legend_title="Gene Category",
    legend = dict(font=dict(size=12), orientation='h', yanchor='top', y=-0.2),
    margin=dict(t=50, b=80))
    
    fig2 = go.Figure()
    fig2.update_layout(
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'))
    fig2.update_layout(plot_bgcolor='white')
    fig2.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5)

    ymax = 0
    for i, cat in enumerate(categories):
        n_hist_graph = go.Histogram(x=df_pred[df_pred['gene_category']==cat][experiment + ' (Z-Score)'].values, name=cat, marker_color=all_colors[color_dict[cat]], opacity=0.50)
        fig2.add_trace(n_hist_graph)
        bins = fig2.full_figure_for_development().data[0].xbins
        nbins = round((bins['end'] - bins['start'])/bins['size'])

        counts, edges = np.histogram(df_pred[df_pred['gene_category']==cat][experiment + ' (Z-Score)'].values, bins=nbins)
        y = max(counts)
        if y > ymax:
            ymax = y
        ymax += 10

    for cat in categories:
        crit_vals = find_alpha_05(df_pred[df_pred['gene_category']==cat][experiment + ' (Z-Score)'].values)
        fig2.add_trace(go.Scatter(x=[crit_vals[0], crit_vals[0], None, crit_vals[1], crit_vals[1]],
                                 y = [0, ymax, None, 0, ymax], 
                                 mode='lines',
                                 line=dict(color=adjust_lightness(all_colors[color_dict[cat]], 0.7), width=2, dash='dot'),
                                 name=cat + ' (alpha=0.05)'),)
    fig2.update_layout(barmode='overlay')
    fig2.update_layout(
    xaxis_title="Predicted Z-Score",
    yaxis_title="Number of Genes",
    legend_title="Gene Category",
    legend = dict(font=dict(size=12), orientation='h', yanchor='top', y=-0.2),
    margin=dict(t=50, b=80))

    fig3 = ff.create_distplot([df_pred[df_pred['gene_category']==cat][experiment + ' (CERES Pred)'].values for cat in categories],
                              categories,
                              colors=[all_colors[color_dict[cat]] for cat in categories],
                              bin_size = 0.1,
                              show_curve = True,
                              show_rug = True)
    fig3.update_layout(
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'))
    fig3.update_layout(plot_bgcolor='white')
    fig3.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5)
    ymax = 0
    fig3.update_layout(barmode='overlay')
    fig3.update_layout(
    xaxis_title="Predicted CERES Score",
    yaxis_title="Probability Density",
    legend_title="Gene Category",
    legend = dict(font=dict(size=12), orientation='h', yanchor='top', y=-0.2),
    margin=dict(t=50, b=80))

    fig4 = ff.create_distplot([df_pred[df_pred['gene_category']==cat][experiment + ' (Z-Score)'].values for cat in categories],
                              categories,
                              colors=[all_colors[color_dict[cat]] for cat in categories],
                              bin_size = 0.1,
                              show_curve = True,
                              show_rug = True)
    fig4.update_layout(
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'))
    fig4.update_layout(plot_bgcolor='white')
    fig4.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5)
    ymax = 0
    fig4.update_layout(barmode='overlay')
    fig4.update_layout(
    xaxis_title="Predicted Z-Scores",
    yaxis_title="Probability Density",
    legend_title="Gene Category",
    legend = dict(font=dict(size=12), orientation='h', yanchor='top', y=-0.2),
    margin=dict(t=50, b=80))
    

    return fig, fig2, fig3, fig4

@app.callback(Output('hist-gene-table', 'children'),
              Input('total_dist_graph', 'clickData'),
              Input('z_score_dist_graph', 'clickData'),
              State('total_dist_graph', 'figure'),
              State('z_score_dist_graph', 'figure'),
              State('df_pred-store', 'data'),
              State('choose_experiment_2', 'value'),
              prevent_initial_call = True)
def update_gene_table(total_data, z_score_data, total_fig, z_score_fig, pred_data, experiment):
    trigger_id = ctx.triggered_id
    if trigger_id == 'total_dist_graph':
        if total_data is not None and pred_data is not None and experiment is not None and total_fig is not None:
            df_pred = pd.DataFrame(pred_data)
            click_data = total_data
            figure = total_fig

            curve = click_data['points'][0]['curveNumber']
            point_numbers = click_data['points'][0]['pointNumbers']
            trace = figure['data'][curve]['name']

            df_pred = df_pred[df_pred['gene_category']==trace]
            df_pred = df_pred.reset_index()
            df_pred = df_pred[df_pred.index.isin(point_numbers)]
            columns = ['gene', experiment + ' (CERES Pred)', experiment + ' (Z-Score)']
            df_pred = df_pred[columns]
            df_pred['gene'] = df_pred['gene'].apply(lambda x: create_gene_link(x, 'uniprot', 'gene'))
            df_pred = df_pred.rename(columns={'gene': 'Gene', experiment + ' (CERES Pred)': 'CERES Pred', experiment + ' (Z-Score)': 'Z-Score'})
            df_pred = df_pred.round(3)

            return dbc.Table.from_dataframe(df_pred, striped=True, bordered=True, hover=True)
        
    if trigger_id == 'z_score_dist_graph':
        if pred_data is not None and z_score_data is not None and experiment is not None and z_score_fig is not None:
            df_pred = pd.DataFrame(pred_data)
            click_data = z_score_data
            figure = z_score_fig

            curve = click_data['points'][0]['curveNumber']
            point_numbers = click_data['points'][0]['pointNumbers']
            trace = figure['data'][curve]['name']

            df_pred = df_pred[df_pred['gene_category']==trace]
            df_pred = df_pred.reset_index()
            df_pred = df_pred[df_pred.index.isin(point_numbers)]
            columns = ['gene', experiment + ' (CERES Pred)', experiment + ' (Z-Score)']
            df_pred = df_pred[columns]
            df_pred['gene'] = df_pred['gene'].apply(lambda x: create_gene_link(x, 'uniprot', 'gene'))
            df_pred = df_pred.rename(columns={'gene': 'Gene', experiment + ' (CERES Pred)': 'CERES Pred', experiment + ' (Z-Score)': 'Z-Score'})
            df_pred = df_pred.round(3)
            return dbc.Table.from_dataframe(df_pred, striped=True, bordered=True, hover=True)

#-----------------------------------------------------------GENES PAGE CALLBACKS------------------------------------------------------------------------------

@app.callback(Output('pick-community-1', 'options'),
          Output('pick-community-2', 'options'),
          Output('pick-community-3', 'options'),
          Input('dummy_div_genes', 'children'),)
def update_dropdown(aux):

    df_clusters = pd.read_csv('./data/feat_summary_varExp_filtered_class.csv')
    landmark_genes = sorted(list(set(df_clusters['feature'].tolist())))
    landmark_genes.insert(0, 'ALL')

    return landmark_genes, landmark_genes, landmark_genes

@app.callback(Output('gene-list-dropdown_1', 'options'),
              Output('gene-list-dropdown_2', 'options'),
              Input('pick-community-1', 'value'),
              prevent_initial_call=True)
def update_gene_list_1(community):
    if community is None:
        return dash.no_update
    else:
        if community == 'ALL':
            with open(DepMap19Q4_GENES_PICKLE, 'rb') as f:
                gene_list = pickle.load(f)
            return gene_list, gene_list
        else:
            clust = df_clusters[df_clusters['feature'] == community]['target'].tolist()
            
            clust2 = clust.copy()
            for i,j in G.edges(clust):
                clust2.append(j)
            clust2 = list(set(clust2))

            return clust2, clust2

@app.callback(Output('gene-list-dropdown_3', 'options'),
              Output('gene-list-dropdown_4', 'options'),
              Input('pick-community-2', 'value'),
              prevent_initial_call=True)
def update_gene_list_1(community):
    if community is None:
        return dash.no_update
    else:
        if community == 'ALL':
            with open(DepMap19Q4_GENES_PICKLE, 'rb') as f:
                gene_list = pickle.load(f)
            return gene_list, gene_list
        else:
            clust = df_clusters[df_clusters['feature'] == community]['target'].tolist()
            
            clust2 = clust.copy()
            for i,j in G.edges(clust):
                clust2.append(j)
            clust2 = list(set(clust2))

            return clust2, clust2

@app.callback(Output('multi-gene-list-dropdown', 'options'),
              Input('pick-community-3', 'value'),
              prevent_initial_call=True)
def update_gene_list_1(community):
    if community is None:
        return dash.no_update
    else:
        if community == 'ALL':
            with open(DepMap19Q4_GENES_PICKLE, 'rb') as f:
                gene_list = pickle.load(f)
            return gene_list
        else:
            clust = df_clusters[df_clusters['feature'] == community]['target'].tolist()
            
            clust2 = clust.copy()
            for i,j in G.edges(clust):
                clust2.append(j)
            clust2 = list(set(clust2))

            return clust2

@app.callback(Output('single_dist_graph_1','figure'),
              Input('gene-list-dropdown_1', 'value'),
              State('df_pred-store', 'data'),
              State('experiment_labels', 'data'),
              prevent_initial_call = True)
def run_inference(gene2visualize, pred_data, experiments):
    if gene2visualize is None or pred_data is None or experiments is None:
        return dash.no_update
    df_pred = pd.DataFrame(pred_data)
    all_colors = {}
    all_colors.update(mcolors.TABLEAU_COLORS)
    standard_colors = mcolors.TABLEAU_COLORS.copy()
    del standard_colors['tab:green']
    del standard_colors['tab:blue']
    color_dict = dict(zip(experiments, standard_colors))
    fig = go.Figure()
    fig.update_layout(
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'))
    fig.update_layout(plot_bgcolor='white')
    fig.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5)
    y, x = np.histogram(dm_data[gene2visualize].values[~np.isnan(dm_data[gene2visualize].values)], bins=40)
    ymax = myround(y.max()+5)
    fig.add_trace(go.Histogram(x=dm_data[gene2visualize].values[~np.isnan(dm_data[gene2visualize].values)], xbins=dict(start=x[0], end=x[-1], size=x[1]-x[0]), name=gene2visualize, marker_color=all_colors['tab:blue'], opacity=0.50))
    for e in experiments:
        fig.add_trace(go.Scatter(x=[df_pred[df_pred['gene'] == gene2visualize][e + ' (CERES Pred)'].values[0], df_pred[df_pred['gene'] == gene2visualize][e + ' (CERES Pred)'].values[0]],
                                 y=[0, ymax],
                                 mode='lines',
                                 line=dict(color=all_colors[color_dict[e]], width=1.5, dash='dot'),
                                 name=e + ' (' + gene2visualize + ')'))
    fig.update_layout(barmode='overlay')
    fig.update_layout(plot_bgcolor='white')
    fig.update_xaxes(range=[-3, 1.5])
    fig.update_layout(
    title_text = 'Distribution of all CERES Scores for ' + gene2visualize,
    title_x = 0.5,
    xaxis_title="CERES Score of " + gene2visualize,
    yaxis_title="Number of Cell lines",
    legend_title="Legend",
    legend = dict(font=dict(size=12), orientation='h', yanchor='top', y=-0.2),
    margin=dict(t=50, b=80))
    return fig

@app.callback(Output('single_dist_graph_2','figure'),
              Input('gene-list-dropdown_2', 'value'),
              State('df_pred-store', 'data'),
              State('experiment_labels', 'data'),
              prevent_initial_call = True)
def run_inference(gene2visualize, pred_data, experiments):
    if gene2visualize is None or pred_data is None or experiments is None:
        return dash.no_update
    df_pred = pd.DataFrame(pred_data)
    all_colors = {}
    all_colors.update(mcolors.TABLEAU_COLORS)
    standard_colors = mcolors.TABLEAU_COLORS.copy()
    del standard_colors['tab:green']
    del standard_colors['tab:blue']
    color_dict = dict(zip(experiments, standard_colors))
    fig = go.Figure()
    fig.update_layout(
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'))
    fig.update_layout(plot_bgcolor='white')
    fig.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5)
    y, x = np.histogram(dm_data[gene2visualize].values[~np.isnan(dm_data[gene2visualize].values)], bins=40)
    ymax = myround(y.max()+5)
    fig.add_trace(go.Histogram(x=dm_data[gene2visualize].values[~np.isnan(dm_data[gene2visualize].values)], xbins=dict(start=x[0], end=x[-1], size=x[1]-x[0]), name=gene2visualize, marker_color=all_colors['tab:green'], opacity=0.50))
    for e in experiments:
        fig.add_trace(go.Scatter(x=[df_pred[df_pred['gene'] == gene2visualize][e + ' (CERES Pred)'].values[0], df_pred[df_pred['gene'] == gene2visualize][e + ' (CERES Pred)'].values[0]],
                                 y=[0, ymax],
                                 mode='lines',
                                 line=dict(color=all_colors[color_dict[e]], width=1.5, dash='dot'),
                                 name=e + ' (' + gene2visualize + ')'))
    fig.update_layout(barmode='overlay')
    fig.update_layout(plot_bgcolor='white')
    fig.update_xaxes(range=[-3, 1.5])
    fig.update_layout(
    title_text = 'Distribution of all CERES Scores for ' + gene2visualize,
    title_x = 0.5,
    xaxis_title="CERES Score of " + gene2visualize,
    yaxis_title="Number of Cell lines",
    legend_title="Legend",
    legend = dict(font=dict(size=12), orientation='h', yanchor='top', y=-0.2),
    margin=dict(t=50, b=80))
    return fig

@app.callback(Output('depmap_link_1', 'children'),
              Output('uniprot_link_1', 'children'),
              Output('tumor_portal_link_1', 'children'),
              Output('ncbi_link_1', 'children'),
              Input('gene-list-dropdown_1', 'value'),
              )
def update_gene_1_links(gene):
    if gene is None:
        return dash.no_update
    else:
        depmap_link = create_gene_link(gene, 'depmap', 'DepMap: ' + str(gene))
        uniprot_link = create_gene_link(gene, 'uniprot', 'UniProt: ' + str(gene))
        tumor_portal_link = create_gene_link(gene, 'tumor_portal', 'Tumor Portal: ' + str(gene))
        ncbi_link = create_gene_link(gene, 'ncbi', 'NCBI: ' + str(gene))
        return [depmap_link], [uniprot_link], [tumor_portal_link], [ncbi_link]
    
@app.callback(Output('depmap_link_2', 'children'),
              Output('uniprot_link_2', 'children'),
              Output('tumor_portal_link_2', 'children'),
              Output('ncbi_link_2', 'children'),
              Input('gene-list-dropdown_2', 'value'),
              )
def update_gene_1_links(gene):
    if gene is None:
        return dash.no_update
    else:
        depmap_link = create_gene_link(gene, 'depmap', 'DepMap: ' + str(gene))
        uniprot_link = create_gene_link(gene, 'uniprot', 'UniProt: ' + str(gene))
        tumor_portal_link = create_gene_link(gene, 'tumor_portal', 'Tumor Portal: ' + str(gene))
        ncbi_link = create_gene_link(gene, 'ncbi', 'NCBI: ' + str(gene))
        return [depmap_link], [uniprot_link], [tumor_portal_link], [ncbi_link]

@app.callback(Output('pie_chart_graph', 'figure'),
              Input('pie_chart_slider', 'value'),
              State('df_pred-store', 'data'),
              prevent_initial_call = True)
def run_pie_chart_output(threshold, pred_data):
    df_pred = pd.DataFrame(pred_data)
    df_pred_sig_genes = df_pred[df_pred['ceres_pred'] <= threshold]
    df_pred_unsig_genes = df_pred[df_pred['ceres_pred'] > threshold]

    fig = make_subplots(rows=1, cols=2, specs=[[{"type": "pie"}, {"type": "pie"}]])

    labels = df_pred_sig_genes['gene_category'].value_counts().index
    values = df_pred_sig_genes['gene_category'].value_counts().values
    fig.add_trace(trace=go.Pie(labels=labels, values=values, name='significant genes', title='Significant Genes'),row=1,col=1)

    labels = df_pred_unsig_genes['gene_category'].value_counts().index
    values = df_pred_unsig_genes['gene_category'].value_counts().values
    fig.add_trace(trace=go.Pie(labels=labels, values=values, name='insignificant genes', title='Insignificant Genes'),row=1,col=2)
    return fig

@app.callback(Output('z_score_plots_all', 'figure'),
              Input('average_slider', 'value'),
              Input('prediction_slider', 'value'),
              Input('z_score_slider', 'value'),
              State('df_pred-store', 'data'),
              prevent_initial_call = True)
def run_z_score_graphs_all(avg, pred, z_score, pred_data):
    df_pred = pd.DataFrame(pred_data)
    df_pred_top = df_pred[df_pred['ceres_pred'] <= pred]
    df_pred_top = df_pred_top[df_pred_top['avg'] >= avg]
    df_pred_top = df_pred_top[df_pred_top['z_score'] < z_score]
    df_pred_top = df_pred_top.sort_values(by=['z_score'], ascending =False)
    gene2visualize = df_pred_top.iloc[-10:]['gene'].tolist()
    x_pred = df_pred_top.iloc[-10:]['ceres_pred'].tolist()
    x_avg = df_pred_top.iloc[-10:]['avg'].tolist()
    dm_data_sub = dm_data[gene2visualize]
    fig = make_subplots(rows=5, cols=2)
    g=9
    for i in range (1,6):
        for j in range (1,3):
            print(dm_data_sub[gene2visualize[g]])
            print(gene2visualize[g])
            fig.add_trace(go.Histogram(x = dm_data_sub[gene2visualize[g]], name=gene2visualize[g]), row=i, col=j)
            fig.add_shape(go.layout.Shape(type='line', xref='x', x0=x_pred[g], y0=0, x1=x_pred[g], y1=100, line=dict(color="Red", width=2, dash="dashdot",)), row=i, col=j)
            fig.add_shape(go.layout.Shape(type='line', xref='x', x0=x_avg[g], y0=0, x1=x_avg[g], y1=100, line=dict(color="Green", width=2, dash="dashdot",)), row=i, col=j)
            g-=1
    return fig

@app.callback(Output('pairwise_gene_comp','figure'),
              Input('gene-list-dropdown_3', 'value'),
              Input('gene-list-dropdown_4', 'value'),
              Input('select_color_by', 'value'),
              Input('select_category', 'value'),
              State('experiment_labels', 'data'),
              State('df_pred-store', 'data'),
              prevent_initial_call = True)
def run_pairwise_comp(gene_1, gene_2, color_by, category_by, experiments, pred_data):
    if gene_1 is None or gene_2 is None or pred_data is None:
        return dash.no_update
    elif color_by is None and category_by is None:
        df_pred = pd.DataFrame(pred_data)
        dm_data_small = dm_data.dropna()
        dm_data_small = dm_data_small.merge(df_sample_data, how='left', right_on='DepMap_ID', left_index=True)
        color_list = ['#7f00ff', '#ff7f00', '#0000ff', '#ffff00', '#007fff', '#7fff00', '#00ffff', '#00ff00']
        color_dict = dict(zip(experiments, color_list))
        fig = go.Figure()

        fig.add_trace(go.Scatter(x=dm_data_small[gene_1], y=dm_data_small[gene_2], text=dm_data_small['stripped_cell_line_name'], mode='markers',
                                 marker=dict(
                                     color='lightgrey',
                                     size=4,
                                     opacity=0.7,
                                     line=dict(
                                         color='DarkSlateGray',
                                         width=1
                                     )), name='other cell-lines'))
        for e in experiments:
            gene_1_pred = df_pred[df_pred['gene'] == gene_1][e + ' (CERES Pred)'].values
            gene_2_pred = df_pred[df_pred['gene'] == gene_2][e + ' (CERES Pred)'].values
            fig.add_trace(go.Scatter(x=gene_1_pred, y=gene_2_pred, mode='markers',
                                     marker=dict(
                                         color=color_dict[e],
                                         size=7,
                                         line=dict(
                                             color='DarkSlateGray',
                                             width=1
                                         )), name=e + ' (CERES Pred)'))

            fig.add_trace(go.Scatter(x=gene_1_pred, y=gene_2_pred, mode='markers', 
                    marker=dict(
                    color=color_dict[e],
                    size=14,
                        gradient=dict(
                        color=color_dict[e],
                        type="radial",
                    ),
                    opacity = 0.4,
                    line=dict(
                        color='DarkSlateGray',
                        width=0.2
                    )), name=e + ' (CERES Pred)', showlegend=False))
                
        rg=LinearRegression()
        rg_res = rg.fit(np.array(dm_data_small[gene_1]).reshape(-1,1), dm_data_small[gene_2])
        Y = rg_res.predict(np.array(dm_data_small[gene_1]).reshape(-1,1))

        fig.add_trace(go.Scatter(x=dm_data_small[gene_1], y=Y, 
                                 mode='lines',
                                 opacity=0.8,
                                 line=dict(color='red', width=1.5, dash='dot')))

        fig.update_layout(
            title_text = 'Pairwise CERES score comparison between ' + gene_1 + ' and ' + gene_2,
            title_x = 0.5,
            xaxis_title=gene_1 + ' (CERES Score)',
            yaxis_title=gene_2 + ' (CERES Score)',
            legend_title="Legend",
            legend = dict(font=dict(size=12), orientation='h', yanchor='top', y=-0.2),
            )

        fig.update_layout(plot_bgcolor='white')
        fig.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=0.8)
        fig.update_xaxes(zeroline=True, zerolinecolor='black', zerolinewidth=0.8)
        fig.update_xaxes(dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray'),)
        fig.update_yaxes(dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray'))

        return fig

    elif color_by is not None:
        df_pred = pd.DataFrame(pred_data)
        dm_data_small = dm_data.dropna()
        dm_data_small = dm_data_small.merge(df_sample_data, how='left', right_on='DepMap_ID', left_index=True)
        color_list = ['#7f00ff', '#ff7f00', '#0000ff', '#ffff00', '#007fff', '#7fff00', '#00ffff', '#00ff00']
        color_dict = dict(zip(experiments, color_list))

        color_bys = dm_data_small[color_by].unique()
        color_by_list = ['#9edae5', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#c49c94', '#aec7e8', '#ffbb78', '#c5b0d5', '#ff9896', '#ff9896']
        if len(color_by_list) < len(color_bys):
            color_by_dict = dict(zip(color_bys, cycle(color_by_list)))
        else:
            color_by_dict = dict(zip(color_bys, color_by_list))

        fig = go.Figure()

        for cat in color_bys:
            gene_1_preds = dm_data_small[dm_data_small[color_by]==cat][gene_1]
            gene_2_preds = dm_data_small[dm_data_small[color_by]==cat][gene_2]
            texts = dm_data_small[dm_data_small[color_by]==cat]['stripped_cell_line_name']
            fig.add_trace(go.Scatter(x=gene_1_preds, y=gene_2_preds, text=texts, mode='markers',
                                     marker=dict(
                                         color=color_by_dict[cat],
                                         size=10,
                                         opacity=0.7,
                                         line=dict(
                                             color='DarkSlateGray',
                                             width=1
                                         )), name=cat))

            if category_by == cat:
                rg=LinearRegression()
                rg_res = rg.fit(np.array(gene_1_preds).reshape(-1,1), gene_2_preds)
                Y = rg_res.predict(np.array(gene_1_preds).reshape(-1,1))

                fig.add_trace(go.Scatter(x=gene_1_preds, y=Y, 
                                     mode='lines',
                                     line=dict(color=adjust_lightness(color_by_dict[cat], 0.5), width=2.5, dash='dot'),
                                     name='Linear (' + cat + ')'))
        for e in experiments:
            gene_1_pred = df_pred[df_pred['gene'] == gene_1][e + ' (CERES Pred)'].values
            gene_2_pred = df_pred[df_pred['gene'] == gene_2][e + ' (CERES Pred)'].values
            fig.add_trace(go.Scatter(x=gene_1_pred, y=gene_2_pred, mode='markers',
                                     marker=dict(
                                         color=color_dict[e],
                                         size=10,
                                         line=dict(
                                             color='DarkSlateGray',
                                             width=1
                                         )), name=e + ' (CERES Pred)'))

            fig.add_trace(go.Scatter(x=gene_1_pred, y=gene_2_pred, mode='markers', 
                    marker=dict(
                    color=color_dict[e],
                    size=15,
                        gradient=dict(
                        color=color_dict[e],
                        type="radial",
                    ),
                    opacity = 0.4,
                    line=dict(
                        color='DarkSlateGray',
                        width=0.2
                    )), name=e + ' (CERES Pred)', showlegend=False))
                
        rg=LinearRegression()
        rg_res = rg.fit(np.array(dm_data_small[gene_1]).reshape(-1,1), dm_data_small[gene_2])
        Y = rg_res.predict(np.array(dm_data_small[gene_1]).reshape(-1,1))

        fig.add_trace(go.Scatter(x=dm_data_small[gene_1], y=Y, 
                                 mode='lines',
                                 line=dict(color='red', width=2.5, dash='dot'),
                                 name='Linear (all cell lines)'))

        fig.update_layout(
            title_text = 'Pairwise CERES score comparison between ' + gene_1 + ' and ' + gene_2,
            title_x = 0.5,
            xaxis_title=gene_1 + ' (CERES Score)',
            yaxis_title=gene_2 + ' (CERES Score)',
            legend_title="Legend",
            legend = dict(font=dict(size=12), orientation='v', xanchor='right', x=-0.2),
            )

        fig.update_layout(plot_bgcolor='white')
        fig.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=0.8)
        fig.update_xaxes(zeroline=True, zerolinecolor='black', zerolinewidth=0.8)
        fig.update_xaxes(dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray'),)
        fig.update_yaxes(dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray'))

        return fig

@app.callback(Output('select_category', 'options'),
              Input('select_color_by', 'value'),
              prevent_initial_call = True)
def return_options(color_by):
    dm_data_small = dm_data.dropna()
    dm_data_small = dm_data_small.merge(df_sample_data, how='left', right_on='DepMap_ID', left_index=True)
    values = dm_data_small[color_by].unique()
    return [{'label': value, 'value': value} for value in values]

@app.callback(Output('multi_gene_comp','figure'),
              Input('multi-gene-list-dropdown', 'value'),
              State('experiment_labels', 'data'),
              State('df_pred-store', 'data'))
def run_multiple_genes(gene_list, experiments, pred_data):
    if gene_list is None or experiments is None or pred_data is None:
        raise PreventUpdate
    elif len(gene_list) < 2 or len(gene_list) > 10:
        raise PreventUpdate
    else:
        n = len(gene_list)
        fig = make_subplots(rows=n, cols=n, shared_yaxes=True, shared_xaxes=True, horizontal_spacing=0.01)
        df_pred = pd.DataFrame(pred_data)
        dm_data_small = dm_data.dropna()
        dm_data_small = dm_data_small.merge(df_sample_data, how='left', right_on='DepMap_ID', left_index=True)
    
        color_list = ['#7f00ff', '#ff7f00', '#0000ff', '#ffff00', '#007fff', '#7fff00', '#00ffff', '#00ff00']
        color_dict = dict(zip(experiments, color_list))
    
        for i in range (1,n+1):
            fig.update_yaxes(title_text=gene_list[i-1], row=i, col=1)
            for j in range (1,n+1):
                rg=LinearRegression()
                rg_res = rg.fit(np.array(dm_data_small[gene_list[j-1]]).reshape(-1,1), dm_data_small[gene_list[i-1]])
                Y = rg_res.predict(np.array(dm_data_small[gene_list[j-1]]).reshape(-1,1))
            
                fig.add_trace(go.Scatter(x=dm_data_small[gene_list[j-1]], 
                                         y=dm_data_small[gene_list[i-1]],
                                         text=dm_data_small['stripped_cell_line_name'],
                                         mode='markers',
                marker=dict(
                size=2,
                color = 'lightgrey',
                opacity = 1,line=dict(
                    color='DarkSlateGray',
                    width=0.2
                ))
                ), row=i, col=j)

                for e in experiments:
            
                    gene_1_pred = df_pred[df_pred['gene'] == [gene_list[j-1]][0]][e + ' (CERES Pred)'].values
                    gene_2_pred = df_pred[df_pred['gene'] == [gene_list[i-1]][0]][e + ' (CERES Pred)'].values
            
                    fig.add_trace(go.Scatter(x=gene_1_pred, y=gene_2_pred, mode='markers', 
                    marker=dict(
                    color=color_dict[e],
                    size=10,
                        gradient=dict(
                        color=color_dict[e],
                        type="radial",
                    ),
                    opacity = 0.4,
                    line=dict(
                        color='DarkSlateGray',
                        width=0.2
                    )), name=e), row=i, col=j)
            
                    fig.add_trace(go.Scatter(x=gene_1_pred, y=gene_2_pred, mode='markers', 
                    marker=dict(
                    color=color_dict[e],
                    size=4,
                    opacity = 1,
                    line=dict(
                        color='DarkSlateGray',
                        width=0.2
                    )), name=e), row=i, col=j)
            
                fig.add_trace(go.Scatter(x=dm_data_small[gene_list[j-1]], y=Y, 
                                         mode='lines',
                                         opacity=0.8,
                                         line=dict(color='red', width=1.5, dash='dot')), row=i, col=j)
            
                if i == n:
                    fig.update_xaxes(title_text=gene_list[j-1], row=i, col=j)

        fig.update_layout(plot_bgcolor='white')
        fig.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=0.8)
        fig.update_xaxes(zeroline=True, zerolinecolor='black', zerolinewidth=0.8)
        fig.update_layout(showlegend=False)
        fig.update_xaxes(dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray'),)
        fig.update_yaxes(dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray'))
        fig['layout'].update(height=1000)

        return fig


@app.callback(Output('output_z_score', 'figure'),
              Input('choose_experiment_3', 'value'),
              Input('gene_cat_checklist_2', 'value'),
              State('df_pred-store', 'data'),
              prevent_initial_call = True)
def update_z_score_graph(experiment, categories, pred_data):
    if experiment is None or categories is None or pred_data is None:
        return dash.no_update
    else:
        df_pred = pd.DataFrame(pred_data)
        all_colors = {}
        all_colors.update(mcolors.TABLEAU_COLORS)
        color_dict = dict(zip(["conditional essential", "common nonessential", "common essential"], mcolors.TABLEAU_COLORS))
        fig = go.Figure()
        for cat in categories:
            fig.add_trace(go.Scatter(x=df_pred[df_pred['gene_category']==cat][experiment + ' (CERES Pred)'].values, 
                         y=df_pred[df_pred['gene_category']==cat][experiment + ' (Z-Score)'].values,
            mode='markers',
            marker=dict(
            color=all_colors[color_dict[cat]],
            size=6,
            opacity = 0.5,
            line=dict(
                color='DarkSlateGray',
                width=1
            )),
            customdata=df_pred[df_pred['gene_category']==cat]['gene'].values,
            hovertemplate='gene:%{customdata}',
            name=cat))
        fig.update_layout(xaxis_range=[-2.2,1])
        fig.update_layout(xaxis_title="CERES Score Prediction", yaxis_title="Z-Score Prediction")
        return fig

@app.callback(Output('select_experiment', 'options'),
              Output('select_experiment_cluster', 'options'),
              State('experiment_labels', 'data'),
              Input('dummy_div_genes', 'children'),)
def update_cell_line_dropdown(experiments, aux):
    if experiments is None:
        return [{'label': value, 'value': value} for value in dm_data.index.tolist()], None
    elif experiments is not None:
        return [{'label': value, 'value': value} for value in dm_data.index.tolist()], [{'label': value, 'value': value} for value in experiments]

@app.callback(Output('cell-line-or-experiment', 'data'),
              Input('select_experiment', 'value'),
              Input('select_experiment_cluster', 'value'),
              prevent_initial_call=True)
def update_choice(cell_line, experiment):
    trigger_id = ctx.triggered_id
    if trigger_id == 'select_experiment':
        return 'cell line'
    elif trigger_id == 'select_experiment_cluster':
        return 'experiment'

@app.callback(Output('network_graph_1', 'figure'),
              Output('louvain-gene-table', 'children'),
              Input('select_community', 'value'),
              Input('cell-line-or-experiment', 'data'),
              State('select_experiment', 'value'),
              State('select_experiment_cluster', 'value'),
              State('df_pred-store', 'data'),
              Input('select_network_layout', 'value'),
              Input('sig_genes_threshold', 'value'),
              Input('z-score_threshold', 'value'),
              prevent_initial_call = True)
def create_network_graph_1(feature, choice, cell_line, experiment, pred_data, layout, c_threshold, z_threshold):
    if feature is None or (cell_line is None and experiment is None):
        return dash.no_update
    else:
        clust = df_clusters[df_clusters['feature'] == feature]['target'].tolist()
        
        clust2 = clust.copy()
        for i,j in G.edges(clust):
            clust2.append(j)
        clust2 = list(set(clust2))

        subgraph = nx.Graph()
        subgraph.add_nodes_from(clust2)
        subgraph.add_edges_from(G.edges(clust))

        if choice == 'cell line':

            df_temp = dm_data.copy()
            df_temp.columns = [i.split(' ',1)[0] for i in df_temp.columns.tolist()]
            df_temp = df_temp[[i for i in clust2 if i in df_temp.columns]]
            df_temp = df_temp.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
            df_temp_T = df_temp.T

            df_return = df_temp_T[cell_line].reset_index().rename(columns={'index':'Gene'})
            stats_data = pd.read_csv(gene_stats_dir, index_col=0)

            z_score_list = []

            for i, gene in enumerate(df_return['Gene'].values.tolist()):
                
                z_score = calc_z_score(stats_data, gene, df_return[cell_line][i])
                z_score_list.append(z_score)
            
            df_return['Z-Score'] = z_score_list

            df_return = df_return.rename(columns={cell_line: 'CERES Score'})
            df_return = df_return.round(3)

            sig_genes = df_return[(df_return['CERES Score'] < c_threshold) & (abs(df_return['Z-Score']) >= z_threshold)]['Gene'].values.tolist()
            logging.warning(sig_genes)
            subtitle = cell_line

            df_return['Gene'] = df_return['Gene'].apply(lambda x: create_gene_link(x, 'uniprot', 'gene'))

        elif choice == 'experiment':

            df_pred = pd.DataFrame(pred_data)
            df_pred = df_pred[['gene', 'gene_category', 'avg', 'std', experiment + ' (CERES Pred)', experiment + ' (Z-Score)']]
            df_pred = df_pred[df_pred['gene'].isin(clust2)]

            df_return = df_pred[['gene', experiment + ' (CERES Pred)', experiment + ' (Z-Score)']]
            
            df_return = df_return.rename(columns={'gene': 'Gene', experiment + ' (CERES Pred)': 'CERES Pred', experiment + ' (Z-Score)': 'Z-Score'})
            df_return = df_return.round(3)

            sig_genes = df_return[(df_return['CERES Pred'] < c_threshold) & (abs(df_return['Z-Score']) >= z_threshold)]['Gene'].values.tolist()
            logging.warning(sig_genes)
            subtitle = experiment + ' (CERES Prediction)'

            df_return['Gene'] = df_return['Gene'].apply(lambda x: create_gene_link(x, 'uniprot', 'gene'))

        subgraph_2 = subgraph.subgraph(sig_genes)

        subgraph_3 = nx.Graph()
        subgraph_3.add_nodes_from(clust2)
        subgraph_3.add_edges_from(subgraph_2.edges(sig_genes))

        pos = nx.spring_layout(subgraph)

        if layout == 'spring':
            pos = nx.spring_layout(subgraph)
        elif layout == 'circular':
            pos = nx.circular_layout(subgraph)
        elif layout == 'random':
            pos = nx.random_layout(subgraph)
        elif layout == 'spectral':
            pos =  nx.spectral_layout(subgraph)
        elif layout == 'fruchterman_reingold':
            pos = nx.fruchterman_reingold_layout(subgraph_3)

        fig=go.Figure()

        Xv=[pos[k][0] for k in clust2]
        Yv=[pos[k][1] for k in clust2]
        Xed=[]
        Yed=[]

        for edge in subgraph.edges(clust2):
            Xed+=[pos[edge[0]][0],pos[edge[1]][0], None]
            Yed+=[pos[edge[0]][1],pos[edge[1]][1], None]

        trace1=go.Scatter(x=Xed,
                    y=Yed,
                    mode='lines',
                    line=dict(color='gray', width=0.5),
                    hoverinfo='none',
                    opacity=0.5,
                    name='Background Edges'
                    )
        trace2=go.Scatter(x=Xv,
                    y=Yv,
                    mode='markers',
                    name='Insignificant Gene Effect',
                    marker=dict(symbol='circle',
                                    size=9,
                                    color='lightblue',
                                    line=dict(color='rgb(50,50,50)', width=0.5)
                                    ),
                    text=clust2,
                    hoverinfo='text'
                    )

        fig.add_trace(trace1)
        fig.add_trace(trace2)

        Xv=[pos[k][0] for k in sig_genes]
        Yv=[pos[k][1] for k in sig_genes]
        Xed=[]
        Yed=[]

        for edge in subgraph_3.edges(sig_genes):
            Xed+=[pos[edge[0]][0],pos[edge[1]][0], None]
            Yed+=[pos[edge[0]][1],pos[edge[1]][1], None]

        trace3=go.Scatter(x=Xed,
                    y=Yed,
                    mode='lines',
                    line=dict(color='red', 
                            width=0.5),
                    hoverinfo='none',
                    name='Significant Edges'
                    )
        trace4=go.Scatter(x=Xv,
                    y=Yv,
                    mode='markers',
                    name='Significant Gene Effect',
                    marker=dict(symbol='circle',
                                    size=12,
                                    color='purple',
                                    line=dict(color='rgb(50,50,50)', width=0.5)
                                    ),
                    text=sig_genes,
                    hoverinfo='text'
                    )

        fig.add_trace(trace3)
        fig.add_trace(trace4)

        fig.update_layout(
            title_text = 'Network Topology of Louvain Community with Landmark: ' + feature + ' <br>' + subtitle,
            title_x=0.5,
            legend_title="Legend",
            legend = dict(font=dict(size=12), orientation='h', yanchor='top', y=-0.2),
            plot_bgcolor='white',
            showlegend=True
            )

        fig.update_yaxes(dict(zeroline=False, showgrid=False, showticklabels=False))
        fig.update_xaxes(dict(zeroline=False, showgrid=False, showticklabels=False))
        fig.update_xaxes(showticklabels=False)
        #fig['layout'].update(height=900, width=900)

        return fig, dbc.Table.from_dataframe(df_return, striped=True, bordered=True, hover=True)

@app.callback(Output('cluster_graph_1', 'figure'),
              Output('select_cluster_2', 'options'),
              Output('select_clust_feats_1', 'options'),
              Output('select_clust_feats_2', 'options'),
              Output('select_clust_feats_3', 'options'),
              Output('select_clust_feats_4', 'options'),
              Input('select_community_2', 'value'),
              Input('select_cluster_layout', 'value'))
def create_cluster_graph(landmark, layout):
    if landmark is None or layout is None:
        return dash.no_update
    elif landmark is not None and layout is not None:

        phate_op_key = "cluster_models/" + landmark.replace(' ', '_') + "_phate_op.pkl"
        phate_op = pickle.load(return_s3_object(s3, s3_bucket_name, phate_op_key))
        
        coords_key = "cluster_models/" + landmark.replace(' ', '_') + "_coords.npy"
        coords = np.load(return_s3_object(s3, s3_bucket_name, coords_key))

        clusters = phate.cluster.kmeans(phate_op, k=7)
        
        cluster_cell_lines_key = "cluster_models/" + landmark.replace(' ', '_') + '_cell_lines.pkl'
        cluster_cell_lines = pickle.load(return_s3_object(s3, s3_bucket_name, cluster_cell_lines_key))

        Y_phate_tsne = coords[int(layout[0])]

        fig = go.Figure()

        fig.update_layout(plot_bgcolor='white')
        fig.update_yaxes(zeroline=True, zerolinecolor='rgb(211,211,211,0.8)', zerolinewidth=0.7, showgrid=False, title=layout[1:] + "2", showticklabels=True)
        fig.update_xaxes(zeroline=True, zerolinecolor='rgb(211,211,211,0.8)', zerolinewidth=0.7, showgrid=False, title=layout[1:] + "1", showticklabels=True)

        unique = set(list(clusters))
        color_by_list = ['#9edae5', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#c49c94', '#aec7e8', '#ffbb78', '#c5b0d5', '#ff9896', '#ff9896']
        color_bys = unique
        if len(color_by_list) < len(color_bys):
            color_by_dict = dict(zip(color_bys, cycle(color_by_list)))
        else:
            color_by_dict = dict(zip(color_bys, color_by_list))

        for clust in unique:
            fig.add_trace(go.Scatter(x=[row[0] for row in Y_phate_tsne[clusters==clust]], 
                                    y=[row[1] for row in Y_phate_tsne[clusters==clust]],
                                    mode='markers',
                                    marker=dict(
                                        color=color_by_dict[clust],
                                        size=8,
                                        opacity=0.7,
                                        line=dict(
                                        color='DarkSlateGray',
                                        width=1
                                        )),
                                    text=[cluster_cell_lines[i] for i in range(len(cluster_cell_lines)-1) if [clusters==clust][0][i]==True],
                                    hoverinfo='text',
                                    name = "cluster " + str(clust)))
        
        fig.update_layout(
        title_text = 'PHATE Clusters for Lovaiun Community: ' + landmark + '<br>' + layout[1:] + ' Projection',
        title_x = 0.5,
        legend_title="Legend",
        legend = dict(font=dict(size=12), orientation='h', yanchor='top', y=-0.2),
        margin=dict(t=80, b=80))

        return fig, [{'label': "cluster " + str(value), 'value': value} for value in unique], [{'label': "cluster " + str(value), 'value': value} for value in unique], [{'label': "cluster " + str(value), 'value': value} for value in unique], [{'label': "cluster " + str(value), 'value': value} for value in unique], [{'label': "cluster " + str(value), 'value': value} for value in unique]

@app.callback(Output('select_cell_cluster_2', 'options'),
              Output('select_experiment_cluster_2', 'options'),
              State('experiment_labels', 'data'),
              Input('dummy_div_genes', 'children'),)
def update_cell_line_dropdown(experiments, aux):
    if experiments is None:
        return [{'label': value, 'value': value} for value in dm_data.index.tolist()], None
    elif experiments is not None:
        return [{'label': value, 'value': value} for value in dm_data.index.tolist()], [{'label': value, 'value': value} for value in experiments]

@app.callback(Output('cell-line-or-experiment-2', 'data'),
              Input('select_cell_cluster_2', 'value'),
              Input('select_experiment_cluster_2', 'value'),
              prevent_initial_call=True)
def update_choice(cell_line, experiment):
    trigger_id = ctx.triggered_id
    if trigger_id == 'select_cell_cluster_2':
        return 'cell line'
    elif trigger_id == 'select_experiment_cluster_2':
        return 'experiment'

@app.callback(Output('cluster-gene-table', 'children'),
              Input('select_cluster_2', 'value'),
              Input('cell-line-or-experiment-2', 'data'),
              Input('select_community_2', 'value'),
              State('select_cell_cluster_2', 'value'),
              State('select_experiment_cluster_2', 'value'),
              State('df_pred-store', 'data'),
              prevent_initial_call = True)
def create_network_graph_1(clust, choice, landmark, cell_line, experiment, pred_data):
    if clust is None or (cell_line is None and experiment is None) or landmark is None:
        return dash.no_update
    else:
        feature_key = "cluster_models/" + landmark.replace(' ', '_') + "_feat_importance_rf.pkl"
        feature_dict = pickle.load(return_s3_object(s3, s3_bucket_name, feature_key))
        importance_dict = feature_dict[int(clust)]
        feat_labels = list(importance_dict.keys())
        importances = list(importance_dict.values())
        importances, feat_labels = zip(*sorted(zip(importances, feat_labels), reverse=True))

        cluster = df_clusters[df_clusters['feature'] == landmark]['target'].tolist()
        
        clust2 = cluster.copy()
        for i,j in G.edges(cluster):
            clust2.append(j)
        clust2 = list(set(clust2))

        feats_list = [i for i in feat_labels if i in clust2]

        if choice == 'cell line':

            df_temp = dm_data.copy()
            df_temp.columns = [i.split(' ',1)[0] for i in df_temp.columns.tolist()]
            df_temp = df_temp[[i for i in feats_list if i in df_temp.columns]]
            df_temp = df_temp.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
            df_temp_T = df_temp.T

            df_return = df_temp_T[cell_line].reset_index().rename(columns={'index':'Gene'})
            stats_data = pd.read_csv(gene_stats_dir, index_col=0)

            z_score_list = []

            for i, gene in enumerate(df_return['Gene'].values.tolist()):
                
                z_score = calc_z_score(stats_data, gene, df_return[cell_line][i])
                z_score_list.append(z_score)
            
            df_return['Z-Score'] = z_score_list

            df_return = df_return.rename(columns={cell_line: 'CERES Score'})
            df_return = df_return.round(3)

            df_return['Gene'] = df_return['Gene'].apply(lambda x: create_gene_link(x, 'uniprot', 'gene'))

        elif choice == 'experiment':

            df_pred = pd.DataFrame(pred_data)
            df_pred = df_pred[['gene', 'gene_category', 'avg', 'std', experiment + ' (CERES Pred)', experiment + ' (Z-Score)']]
            df_pred = df_pred[df_pred['gene'].isin(feats_list)]

            df_return = df_pred[['gene', experiment + ' (CERES Pred)', experiment + ' (Z-Score)']]
            
            df_return = df_return.rename(columns={'gene': 'Gene', experiment + ' (CERES Pred)': 'CERES Pred', experiment + ' (Z-Score)': 'Z-Score'})
            df_return = df_return.round(3)

            df_return['Gene'] = df_return['Gene'].apply(lambda x: create_gene_link(x, 'uniprot', 'gene'))
        
        return dbc.Table.from_dataframe(df_return, striped=True, bordered=True, hover=True)

@app.callback(Output('cluster_feats_1', 'figure'),
              Input('select_clust_feats_1', 'value'),
              State('select_community_2', 'value'))
def update_clust_feats_1(clust, landmark):
    if landmark is None and clust is None:
        return dash.no_update
    elif landmark is not None and clust is not None:

        feature_key = "cluster_models/" + landmark.replace(' ', '_') + "_feat_importance_rf.pkl"
        feature_dict = pickle.load(return_s3_object(s3, s3_bucket_name, feature_key))
        importance_dict = feature_dict[int(clust)]
        feat_labels = list(importance_dict.keys())
        importances = list(importance_dict.values())
        importances, feat_labels = zip(*sorted(zip(importances, feat_labels), reverse=True))

        fig = go.Figure()

        fig.update_layout(plot_bgcolor='white')
        fig.update_yaxes(zeroline=True, zerolinecolor='rgb(211,211,211,0.8)', zerolinewidth=0.7, showgrid=False, showticklabels=True)
        fig.update_xaxes(zeroline=True, zerolinecolor='rgb(211,211,211,0.8)', zerolinewidth=0.7, showgrid=False, showticklabels=True)

        unique = set(list(['0', '1', '2', '3', '4', '5', '6']))
        color_by_list = ['#9edae5', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#c49c94', '#aec7e8', '#ffbb78', '#c5b0d5', '#ff9896', '#ff9896']
        color_bys = unique
        if len(color_by_list) < len(color_bys):
            color_by_dict = dict(zip(color_bys, cycle(color_by_list)))
        else:
            color_by_dict = dict(zip(color_bys, color_by_list))

        fig.add_trace(go.Bar(
            y=feat_labels[0:20],
            x=importances[0:20],
            name='Top Features for Cluster' + clust,
            orientation='h',
            marker=dict(
                color=color_by_dict[clust],
                opacity = 0.7,
                line=dict(
                    color='DarkSlateGray',
                    width=1
                    ),)
            ))

        fig.update_layout(
            xaxis_title='Feature Score',
            yaxis_title='Feature Name',
            title_text='Top 20 Features for Cluster ' + str(clust),
            title_x=0.5)

    return fig

@app.callback(Output('cluster_feats_2', 'figure'),
              Input('select_clust_feats_2', 'value'),
              State('select_community_2', 'value'))
def update_clust_feats_1(clust, landmark):
    if landmark is None and clust is None:
        return dash.no_update
    elif landmark is not None and clust is not None:
        feature_key = "cluster_models/" + landmark.replace(' ', '_') + "_feat_importance_rf.pkl"
        feature_dict = pickle.load(return_s3_object(s3, s3_bucket_name, feature_key))
        importance_dict = feature_dict[int(clust)]
        feat_labels = list(importance_dict.keys())
        importances = list(importance_dict.values())
        importances, feat_labels = zip(*sorted(zip(importances, feat_labels), reverse=True))

        fig = go.Figure()

        fig.update_layout(plot_bgcolor='white')
        fig.update_yaxes(zeroline=True, zerolinecolor='rgb(211,211,211,0.8)', zerolinewidth=0.7, showgrid=False, showticklabels=True)
        fig.update_xaxes(zeroline=True, zerolinecolor='rgb(211,211,211,0.8)', zerolinewidth=0.7, showgrid=False, showticklabels=True)

        unique = set(list(['0', '1', '2', '3', '4', '5', '6']))
        color_by_list = ['#9edae5', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#c49c94', '#aec7e8', '#ffbb78', '#c5b0d5', '#ff9896', '#ff9896']
        color_bys = unique
        if len(color_by_list) < len(color_bys):
            color_by_dict = dict(zip(color_bys, cycle(color_by_list)))
        else:
            color_by_dict = dict(zip(color_bys, color_by_list))
            
        fig.add_trace(go.Bar(
            y=feat_labels[0:20],
            x=importances[0:20],
            name='Top Features for Cluster' + clust,
            orientation='h',
            marker=dict(
                color=color_by_dict[clust],
                opacity = 0.7,
                line=dict(
                    color='DarkSlateGray',
                    width=1
                    ),)
            ))
        
        fig.update_layout(
            xaxis_title='Feature Score',
            yaxis_title='Feature Name',
            title_text='Top 20 Features for Cluster ' + str(clust),
            title_x=0.5)

    return fig

@app.callback(Output('cluster_feats_3', 'figure'),
              Input('select_clust_feats_3', 'value'),
              State('select_community_2', 'value'))
def update_clust_feats_1(clust, landmark):
    if landmark is None and clust is None:
        return dash.no_update
    elif landmark is not None and clust is not None:
        feature_key = "cluster_models/" + landmark.replace(' ', '_') + "_feat_importance_rf.pkl"
        feature_dict = pickle.load(return_s3_object(s3, s3_bucket_name, feature_key))
        importance_dict = feature_dict[int(clust)]
        feat_labels = list(importance_dict.keys())
        importances = list(importance_dict.values())
        importances, feat_labels = zip(*sorted(zip(importances, feat_labels), reverse=True))

        fig = go.Figure()

        fig.update_layout(plot_bgcolor='white')
        fig.update_yaxes(zeroline=True, zerolinecolor='rgb(211,211,211,0.8)', zerolinewidth=0.7, showgrid=False, showticklabels=True)
        fig.update_xaxes(zeroline=True, zerolinecolor='rgb(211,211,211,0.8)', zerolinewidth=0.7, showgrid=False, showticklabels=True)

        unique = set(list(['0', '1', '2', '3', '4', '5', '6']))
        color_by_list = ['#9edae5', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#c49c94', '#aec7e8', '#ffbb78', '#c5b0d5', '#ff9896', '#ff9896']
        color_bys = unique
        if len(color_by_list) < len(color_bys):
            color_by_dict = dict(zip(color_bys, cycle(color_by_list)))
        else:
            color_by_dict = dict(zip(color_bys, color_by_list))
            
        fig.add_trace(go.Bar(
            y=feat_labels[0:20],
            x=importances[0:20],
            name='Top Features for Cluster' + clust,
            orientation='h',
            marker=dict(
                color=color_by_dict[clust],
                opacity = 0.7,
                line=dict(
                    color='DarkSlateGray',
                    width=1
                    ),)
            ))

        fig.update_layout(
            xaxis_title='Feature Score',
            yaxis_title='Feature Name',
            title_text='Top 20 Features for Cluster ' + str(clust),
            title_x=0.5)

    return fig

@app.callback(Output('cluster_feats_4', 'figure'),
              Input('select_clust_feats_4', 'value'),
              State('select_community_2', 'value'))
def update_clust_feats_1(clust, landmark):
    if landmark is None and clust is None:
        return dash.no_update
    elif landmark is not None and clust is not None:
        feature_key = "cluster_models/" + landmark.replace(' ', '_') + "_feat_importance_rf.pkl"
        feature_dict = pickle.load(return_s3_object(s3, s3_bucket_name, feature_key))
        importance_dict = feature_dict[int(clust)]
        feat_labels = list(importance_dict.keys())
        importances = list(importance_dict.values())
        importances, feat_labels = zip(*sorted(zip(importances, feat_labels), reverse=True))

        fig = go.Figure()

        fig.update_layout(plot_bgcolor='white')
        fig.update_yaxes(zeroline=True, zerolinecolor='rgb(211,211,211,0.8)', zerolinewidth=0.7, showgrid=False, showticklabels=True)
        fig.update_xaxes(zeroline=True, zerolinecolor='rgb(211,211,211,0.8)', zerolinewidth=0.7, showgrid=False, showticklabels=True)

        unique = set(list(['0', '1', '2', '3', '4', '5', '6']))
        color_by_list = ['#9edae5', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#c49c94', '#aec7e8', '#ffbb78', '#c5b0d5', '#ff9896', '#ff9896']
        color_bys = unique
        if len(color_by_list) < len(color_bys):
            color_by_dict = dict(zip(color_bys, cycle(color_by_list)))
        else:
            color_by_dict = dict(zip(color_bys, color_by_list))
            
        fig.add_trace(go.Bar(
            y=feat_labels[0:20],
            x=importances[0:20],
            name='Top Features for Cluster' + clust,
            orientation='h',
            marker=dict(
                color=color_by_dict[clust],
                opacity = 0.7,
                line=dict(
                    color='DarkSlateGray',
                    width=1
                    ),)
            ))
                
        fig.update_layout(
            xaxis_title='Feature Score',
            yaxis_title='Feature Name',
            title_text='Top 20 Features for Cluster ' + str(clust),
            title_x=0.5)

    return fig

app.register_celery_tasks()

if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8080, debug=True)
