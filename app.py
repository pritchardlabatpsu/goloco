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
from scipy.stats import pearsonr, norm
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
from dash.dash_table.Format import Format, Scheme
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from app_session import infer
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.linear_model import LinearRegression
import dash_cytoscape as cyto
import math
import requests
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
from bisect import bisect
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
cell_line_DIRECTORY = './data/cell_line_names.pkl'
manual_DIRECTORY = './data/goloco_user_manual.pdf'

launch_uid = uuid.uuid4()

# AWS Bucket Information
#aws_available = os.environ['AWS_S3_AVAILABLE']

#s3_client = boto3.client('s3')
#s3_bucket_name = 'ceres-infer'
#s3 = boto3.resource('s3', aws_access_key_id=os.environ['AWS_ACCESS_KEY_ID'],
#                          aws_secret_access_key=os.environ['AWS_SECRET_ACCESS_KEY'])
#my_bucket = s3.Bucket(s3_bucket_name)

# read depmap data from pickled directories
global dm_data
dm_data = pd.read_pickle(DepMap19Q4_DIRECTORY)
df_sample_data = pd.read_pickle(DepMap19Q4_SAMPLE_INFO)

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

# define dash app
app = DashProxy(__name__, use_pages=True, suppress_callback_exceptions=True,
                background_callback_manager=background_callback_manager,
                external_stylesheets=[dbc.themes.BOOTSTRAP],
                meta_tags=[{'name': 'viewport', 'content': 'width=device-width, initial-scale=1'}], 
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
                dbc.DropdownMenuItem("Hits", href="/overview"),
                dbc.DropdownMenuItem("Regressions", href="/genes"),
            ],
            nav=True,
            in_navbar=True,
            label="Explore",
        ),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("gProfiler", href="/enrichment"),
                dbc.DropdownMenuItem("Communities", href="/clusters"),
            ],
            nav=True,
            in_navbar=True,
            label="Pathways",
        ),
    ],
    brand="goloco",
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
    dcc.Store(id='z-score-hits', storage_type='memory'),
    dcc.Store(id='p-value-hits', storage_type='memory'),
    dcc.Store(id='enrichment-hits', storage_type='session'),
    dcc.Loading(id='store-data-loading', children = [html.Div([dcc.Store(id='df_pred-store', storage_type='session')])], type='circle', fullscreen=True),
    dcc.Store(id='df_pred-store-tmp-1', storage_type='memory'),
    dcc.Store(id='df_pred-store-tmp-2', storage_type='memory'),
    dcc.Store(id='df_pred-store-tmp-3', storage_type='memory'),
    dcc.Store(id='experiment_labels', storage_type='session'),
    dcc.Store(id='df-counts-convert-tmp', storage_type='memory'),
    dcc.Store(id='df-counts-convert', storage_type='session'),
    dcc.Store(id='cell-line-or-experiment', storage_type='memory'),
    dcc.Store(id='cell-line-or-experiment-2', storage_type='memory'),
    dcc.Location(id='url_name'),
    dcc.Download(id='download-pred'),
    dcc.Download(id='download-convert'),
    dcc.Download(id='download-z-genes'),
    dcc.Download(id='download-p-genes'),

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
            dbc.ModalHeader("Completed"),
            dbc.ModalBody("Your previous prediction has been uploaded. Select the button below to continue. Navigate to the Explore tabs to visualize your results."),
            dbc.ModalFooter(
                dbc.Button(
                    "Continue", id="completed-tmp-2", color="primary", className="ms-auto", n_clicks=0
                )
            ),
        ],
        id="completed-2-popup",
        is_open=False,
    ),

    dbc.Modal(
        [
            dbc.ModalHeader("Completed"),
            dbc.ModalBody("Your genome wide inference has successfully completed. Select the button below to continue. Navigate to the Explore tabs to visualize the results."),
            dbc.ModalFooter(
                dbc.Button(
                    "Continue", id="completed-tmp-1", color="primary", className="ms-auto", n_clicks=0
                )
            ),
        ],
        id="completed-1-popup",
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

    dbc.Modal(
        [
            dbc.ModalBody("Select 'continue' to run gProfiler functional enrichment analysis and you will be redirected to another page."),
            dbc.ModalFooter(
                dbc.Row([dbc.Col(dbc.Button(
                    "Cancel", id="stay-overview", color="primary", className="ms-auto", n_clicks=0)), 
                         dbc.Col(dbc.Button(
                    "Continue", id="continue-gprofile", color="success", className="ms-auto", n_clicks=0
                ))]),
            ),
        ],
        id="gprofile-overview-popup",
        is_open=False,
    ),

    dcc.Location(id='url', refresh=True),

    html.Div(id='dummy_div_app'), 
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
    elif website == 'ensemble':
        url = generate_ensemble_url(gene)

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

def generate_ensemble_url(gene):
    return f'https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g={gene}'

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

def logScale(input, input_start=1, input_end=50000, output_start=2, output_end=10):
    m = (output_end - output_start) / (np.log(input_end) - np.log(input_start))
    b = -m * np.log(input_start) + output_start
    output = m * np.log(input) + b
    return output

def xScale(input, input_start=1, input_end=1000, output_start=2, output_end=200):
    m = (output_end - output_start) / (input_end - input_start)
    b = -m * input_start + output_start
    output = m * input + b
    return output

def summarize_meta(meta):
    dcol1 = ['version', 'date', 'organism', 'all results', 'ordered', 'no iea', 'sources', 'multiquery', 
              'numeric ns', 'domain scope', 'measure underrepresentation', 'significance threshold method', 
              'user threshold', 'no evidences', 'highlight results']
    dcol2 = [meta['version'], meta['timestamp'], meta['query_metadata']['organism'],
              str(meta['query_metadata']['all_results']),
              str(meta['query_metadata']['ordered']),
              str(meta['query_metadata']['no_iea']),
              ', '.join([str(i) for i in meta['query_metadata']['sources']]),
              str(meta['query_metadata']['combined']),
              meta['query_metadata']['numeric_ns'],
              meta['query_metadata']['domain_scope'],
              str(meta['query_metadata']['measure_underrepresentation']),
              meta['query_metadata']['significance_threshold_method'],
              meta['query_metadata']['user_threshold'],
              str(meta['query_metadata']['no_evidences']),
              str(meta['query_metadata']['highlight'])]
    
    d_details=pd.DataFrame([dict(zip(dcol1, dcol2))])

    return d_details

def color(style, text):
    return f"<span style='font-size:{str(style[1])}'><span style='color:{str(style[0])}'> {str(text)} </span></span>"

def gostplot(gostres, query=False, capped=True, pal=None):
    # gostres is the GOSt response list (contains results and metadata)
    # This function will plot only the sources that were asked from the query

    if pal is None:
        pal = {
            "GO:MF": "#dc3912",
            "GO:BP": "#ff9900",
            "GO:CC": "#109618",
            "KEGG": "#dd4477",
            "REAC": "#3366cc",
            "WP": "#0099c6",
            "TF": "#5574a6",
            "MIRNA": "#22aa99",
            "HPA": "#6633cc",
            "CORUM": "#66aa00",
            "HP": "#990099"
        }

    if 'result' not in gostres:
        raise ValueError("Name 'result' not found from the input")
    if 'meta' not in gostres:
        raise ValueError("Name 'meta' not found from the input")

    df = pd.DataFrame(gostres['result'])
    meta = gostres['meta']
    essential_names = ["source_order", "term_size", "term_name", "term_id", "source", "significant"]

    if query:
        chosen_sources = query
        
        widths = [meta['result_metadata'][source]['number_of_terms'] for source in chosen_sources]
        widthscale = {source: meta['result_metadata'][source]['number_of_terms'] for source in chosen_sources}
        
        space = 1000
        starts = []
        starts.append(1000)

        if not len(widthscale) < 2:
            for i in range(1, len(widthscale)):
                starts.append(starts[i-1] + space + widths[i - 1])
                
        starts = dict(zip(chosen_sources, starts))

        sourcediff = set(chosen_sources) - set(pal.keys())
        colors = [color for color in plt.cm.tab20.colors if color not in pal.values()]
        
        df = df[df['source'].isin(chosen_sources)]
    
    else:
        chosen_sources = meta['query_metadata']['sources']
        
        widthscale = {source: meta['result_metadata'][source]['number_of_terms'] for source in chosen_sources}
        starts = {source: sum(widthscale[s] for s in chosen_sources[:i]) + 1000 * i for i, source in enumerate(chosen_sources)}
        
        sourcediff = set(chosen_sources) - set(pal.keys())
        colors = [color for color in plt.cm.tab20.colors if color not in pal.values()]
        
    if sourcediff:
        use_cols = np.random.choice(colors, len(sourcediff), replace=False)
        pal.update(dict(zip(sourcediff, use_cols)))

    df['term_size_scaled'] = df['term_size'].apply(lambda x: logScale(x))
    df['order'] = df.apply(lambda row: xScale(row['source_order'], input_start = 1, input_end=widthscale[row['source']], output_start=starts[row['source']], output_end = starts[row['source']] + widthscale[row['source']]), axis=1)
    df['order'] = df['order'].apply(lambda x: xScale(x))
    df['logpval'] = df['p_value'].apply(lambda x: -math.log10(x))
    
    if capped:
        df['logpval'].loc[df['logpval'] > 16] = 15
        ymin, ymax = -1, 18.5
        ticklabels = ["0", "2", "4", "6", "8", "10", "12", "14", ">16"]
        tickvals = [0, 2, 4, 6, 8, 10, 12, 14, 16]
    else:
        ymax = round(max(df['logpval']))
        ymin = -1

    xticklabels = [source + ' (' + str(df[df['source']==source].shape[0]) + ')' for source in chosen_sources]

    fig = go.Figure()

    for source in chosen_sources:
        xstart = xScale(starts[source])
        xend = xScale(starts[source] + widthscale[source])
        
        xmax = xend

        source_df = df[df['source'] == source]
        fig.add_trace(go.Scatter(
            x=source_df['order'],
            y=source_df['logpval'],
            mode='markers',
            marker=dict(
                color=pal[source],
                size=2*source_df['term_size_scaled'],
                opacity=np.where(source_df['significant'], 0.8, np.where(source_df['p_value'] == 1, 0, 0.2)),
                line=dict(width=0),
            ),
            name=source,
            customdata=source_df[['native', 'term_size', 'name', 'p_value']],
            hovertemplate='<b>%{customdata[0]} (%{customdata[1]})</b><br>%{customdata[2]}<br>%{customdata[3]:.3e}<extra></extra>',
            texttemplate='%{customdata[0]}'
        ))
        
        source_df = df[(df['source'] == source) & (df['highlighted']==True)]
        fig.add_trace(go.Scatter(
            x=source_df['order'],
            y=source_df['logpval'],
            mode='markers',
            marker=dict(
                color=pal[source],
                size=2*source_df['term_size_scaled'],
                opacity=np.where(source_df['significant'], 0.8, np.where(source_df['p_value'] == 1, 0, 0.2)),
                line=dict(width=2,
                          color='black'),
            ),
            name=source,
            customdata=source_df[['native', 'term_size', 'name', 'p_value']],
            hovertemplate='<b>%{customdata[0]} (%{customdata[1]})</b><br>%{customdata[2]}<br>%{customdata[3]:.3e}<extra></extra>',
        ))    

        source_df = source_df[['order', 'logpval', 'native', 'term_size', 'name', 'p_value', 'group_id']]
        '''for row in source_df.iterrows():
                fig.add_annotation(x=row[1]['order'], y=row[1]['logpval'],
                text='<b>({})</b>'.format(row[1]['group_id']),
                showarrow=True,
                font=dict(size=12,
                        color='black'))'''
            
        fig.add_shape(type="line",
                        x0=xstart,
                        x1=xend,
                        y0=-0.5,
                        y1=-0.5,
                        line=dict(color=pal[source], 
                        width=10))

    fig.add_shape(type="line",
                    x0=0,
                    x1=0,
                    y0=0,
                    y1=16,
                    line=dict(color='black', 
                    width=4))

    fig.add_shape(type="line",
                    x0=0,
                    x1=xmax,
                    y0=-0.77,
                    y1=-0.77,
                    line=dict(color='black', 
                    width=2))
    
    '''fig.add_shape(type="line",
                x0=0,
                x1=xmax,
                y0=16,
                y1=16,
                line=dict(color='lightgray', 
                          width=1.5,
                         dash='dot'))'''

    fig.update_layout(
        plot_bgcolor='white',
        showlegend=False,
        title_text = '<b>gProfiler2 Enrichment of Functional Categories',
        title_y = 0.9,
        title_x = 0.55,
        font_family = 'Arial',
        font_color = 'black',
        font_size=16,
        yaxis_title="<b>-log10(p-adj)</b>",
        yaxis=dict(
            range=[ymin, ymax],
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticklabels,
            ticks="outside", 
            tickwidth=1.5,
            ticklen=5,
            tickfont=dict(size=14),
            zeroline=False,
            showgrid=False,
           # title_standoff = 0,
            showspikes=True,
            spikethickness=1,
            tickfont_family="Arial Black",
        ),
        xaxis_title=None,
        xaxis=dict(
            tickvals=[(xScale(starts[source]) + xScale(starts[source] + widthscale[source])) / 2 for source in chosen_sources],
            ticktext=xticklabels,
            tick0 = -1,
            range=[-4, xmax],
            #rangemode = "nonnegative",
            zeroline=False,
            tickangle=-30,
            tickfont=dict(size=14),
            tickfont_family="Arial Black",
        ),
        margin=dict(t=40, r=0, b=20, l=0),
        height=400,
        hoverlabel=dict(
            bgcolor='whitesmoke',
            font_size=12,
            font_family="Arial",
            align="auto",
            ),
    )
    
    fig.add_annotation(x=xmax, y=16,
        text="values above this threshold are capped",
        showarrow=False,
        xshift=-100,
        yshift=-10,
        font=dict(size=10,
                  color='darkgray'))
    
    highlighted_df = df[df['highlighted']==True].sort_values(by='group_id')
    highlighted_df = highlighted_df[['group_id', 'source', 'native', 'name', 'term_size', 'p_value', 'p_value']]
    highlighted_df.columns = ['group_id', 'source', 'native', 'name', 'term_size', 'p_value_2', 'p_value']
    highlighted_df['p_value_2'] = highlighted_df['p_value_2'].apply(lambda x: math.log(x))
    
    more_styles = discrete_background_color_bins(highlighted_df, n_bins=100, columns=['p_value_2'])
    
    table = dash_table.DataTable(data=highlighted_df.to_dict('records'),
                                        columns=[{'name': 'ID', 'id': 'group_id', 'type':'text'},
                                                {'name': 'Source', 'id': 'source', 'type': 'text'},
                                                {'name': 'Term ID', 'id': 'native', 'type': 'text'},
                                                {'name': 'Term Name', 'id': 'name', 'type': 'text'},
                                                {'name': 'Term Size', 'id': 'term_size', 'type': 'text'},
                                                {'name': '', 'id': 'p_value_2', 'type': 'numeric', 'format': Format(precision=2, scheme=Scheme.exponent)},
                                                {'name': 'p-value (adj)', 'id': 'p_value', 'type': 'numeric', 'format': Format(precision=2, scheme=Scheme.exponent)}],
                                        filter_action='native',
                                        sort_action='native',
                                        sort_mode='multi',
                                        style_table={
                                            'height': 600,
                                            'overflowX': 'auto',
                                            'overflowY': 'auto'
                                        },
                                        style_data={
                                            'whiteSpace': 'normal',
                                            'height': 'auto',
                                        },
                                        style_header={
                                                    'backgroundColor': 'lightgrey',
                                                    'font_family': 'Arial',
                                                    'font_size': '16px',
                                                    'paddingLeft': '20px',
                                                    'color': 'black',
                                                    'fontWeight': 'bold'
                                                },
                                        style_cell={'textAlign': 'left',                               
                                                    'font_family': 'Arial',
                                                    'font_size': '16px',
                                                    'paddingLeft': '20px',
                                                    'overflow': 'hidden',
                                                    'textOverflow': 'ellipsis',
                                                    },
                                        style_data_conditional=more_styles,)
    df_meta = summarize_meta(meta)
    df_mapping = pd.DataFrame(meta['genes_metadata']['query']['query_1']['mapping']).T.reset_index().rename(columns={'index': 'HGNC ID', 0: 'Ensembl ID'})
    df_mapping['HGNC ID'] = df_mapping['HGNC ID'].apply(lambda x: create_gene_link(x, 'uniprot', 'gene'))
    df_mapping['Ensembl ID'] = df_mapping['Ensembl ID'].apply(lambda x: create_gene_link(x, 'ensemble', 'gene'))
    return fig, table, [dbc.Table.from_dataframe(df_mapping, striped=True, bordered=True, index=False, size='sm')], [dbc.Table.from_dataframe(df_meta, striped=False, bordered=False, index=False, size='sm')]

def discrete_background_color_bins(df, n_bins=100, columns='all'):
    bounds = [i * (1.0 / n_bins) for i in range(n_bins + 1)]
    if columns == 'all':
        if 'id' in df:
            df_numeric_columns = df.select_dtypes('number').drop(['id'], axis=1)
        else:
            df_numeric_columns = df.select_dtypes('number')
    else:
        df_numeric_columns = df[columns]
    df_max = df_numeric_columns.max().max()
    df_min = df_numeric_columns.min().min()
    ranges = [
        ((df_max - df_min) * i) + df_min
        for i in bounds
    ]
    styles = [{
                                            'if': {'row_index': 'odd'},
                                            'backgroundColor': 'whitesmoke',
                                            },
                                        {'if': {
                                                'filter_query': '{source} = "GO:MF"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#dc3912",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                            
                                          {'if': {
                                                'filter_query': '{source} = "GO:BP"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#ff9900",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                            
                                          {'if': {
                                                'filter_query': '{source} = "GO:CC"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#109618",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                            
                                          {'if': {
                                                'filter_query': '{source} = "KEGG"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#dd4477",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                            
                                          {'if': {
                                                'filter_query': '{source} = "REAC"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#3366cc",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                            
                                          {'if': {
                                                'filter_query': '{source} = "WP"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#0099c6",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                            
                                          {'if': {
                                                'filter_query': '{source} = "TF"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#5574a6",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                            
                                          {'if': {
                                                'filter_query': '{source} = "MIRNA"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#22aa99",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                            
                                          {'if': {
                                                'filter_query': '{source} = "HPA"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#6633cc",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                            
                                          {'if': {
                                                'filter_query': '{source} = "CORUM"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#66aa00",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                            
                                          {'if': {
                                                'filter_query': '{source} = "HP"',
                                                'column_id': ['source', 'native']
                                            },
                                            'backgroundColor': "#990099",
                                            'color': 'white',
                                            'fontWeight': 'bold',},
                                        {'if': {'column_id': 'p_value_2'},
                                                 'width': '1%'},
                            
                                        ]
    colors = sns.color_palette('viridis', n_colors=len(bounds))
    colors = colors.as_hex()
    
    for i in range(1, len(bounds)):
        min_bound = ranges[i - 1]
        max_bound = ranges[i]
        backgroundColor = colors[i - 1]
        color = 'white' if i > len(bounds) / 2. else 'inherit'

        for column in df_numeric_columns:
            styles.append({
                'if': {
                    'filter_query': (
                        '{{{column}}} >= {min_bound}' +
                        (' && {{{column}}} < {max_bound}' if (i < len(bounds) - 1) else '')
                    ).format(column=column, min_bound=min_bound, max_bound=max_bound),
                    'column_id': column
                },
                'backgroundColor': backgroundColor,
                'color': backgroundColor
            })

    return styles

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
              Input('df_pred-store-tmp-3', 'data'),
              Input('completed', 'n_clicks'),
              State('completed-popup', 'is_open'),
              prevent_initial_call = True)
def toggle_modal(n1, n2, is_open):
    trigger_id = ctx.triggered_id
    if n1 or n2:
        return not is_open
    return is_open

# work
@app.callback(Output('completed-1-popup', 'is_open'),
              Input('df_pred-store-tmp-1', 'data'),
              Input('completed-tmp-1', 'n_clicks'),
              State('completed-1-popup', 'is_open'),
              prevent_initial_call = True)
def toggle_modal(n1, n2, is_open):
    trigger_id = ctx.triggered_id
    if n1 or n2:
        return not is_open
    return is_open

@app.callback(Output('completed-2-popup', 'is_open'),
              Input('df_pred-store-tmp-2', 'data'),
              Input('completed-tmp-2', 'n_clicks'),
              State('completed-2-popup', 'is_open'),
              prevent_initial_call = True)
def toggle_modal(n1, n2, is_open):
    trigger_id = ctx.triggered_id
    if n1 or n2:
        return not is_open
    return is_open
# 

@app.callback(Output('conversion-completed-popup', 'is_open'),
              Input('df-counts-convert-tmp', 'data'),
              Input('conversion_completed', 'n_clicks'),
              State('conversion-completed-popup', 'is_open'),
              prevent_initial_call = True)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open

@app.callback(
    Output('gprofile-overview-popup', 'is_open'),
    [Input('stay-overview', 'n_clicks'), 
    Input('continue-gprofile', 'n_clicks'),
    Input('run-gprofiler-zscores', 'n_clicks'),
    Input('run-gprofiler-pvalue', 'n_clicks')],
    [State('gprofile-overview-popup', 'is_open')],
    prevent_initial_call = True
)
def toggle_modal(n1, n2, n3, n4, is_open):
    if n1 or n2 or n3 or n4:
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
            l200_genes = pd.read_pickle(L200_GENES)
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

        # for s3:
        #pkl_filename = "chronos/" + user_id + "_l200_gene_effect_conversion.pkl"
        #pickle_buffer = io.BytesIO()
        #gene_effects.to_pickle(pickle_buffer)
        #s3.Object(s3_bucket_name, pkl_filename).put(Body=pickle_buffer.getvalue())

        # for local:
        pkl_filename = "./ceres-infer/gene_effects/" + user_id + "_l200_gene_effect_conversion.pkl"
        gene_effects.to_pickle(pkl_filename)

        return pkl_filename, None

@app.callback(ServersideOutput('df-counts-convert', 'data'),
              Input('df-counts-convert-tmp', 'data'),
              prevent_initial_call = True) #check
def update_df_convert(pkl_filename):
    if pkl_filename is None:
        return dash.no_update
    else:
        # for s3:
        #obj = s3.Object(s3_bucket_name, pkl_filename)
        #data = obj.get()['Body'].read()
        #gene_effects = pickle.load(io.BytesIO(data))                 
        
        # for local:
        with open(pkl_filename, 'rb') as f:
            gene_effects = pickle.load(f)
        
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
        z_scores = np.zeros(shape=(len(scope),len(a.experiments)))
        x_avgs = []
        x_stds = []
        t1 = time.perf_counter()
        for i, gene in enumerate(scope):
            logging.warning(gene)
            pred = a.infer_gene(gene)
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
        #df_pred = pd.read_csv(os.path.join(OUTPUT_DIRECTORY, 'prediction.csv')) # for testing
        
        # for s3:
        #pkl_filename = "predictions/" + user_id + "_genome_wide_prediction.pkl"
        #pickle_buffer = io.BytesIO()
        #df_pred.to_pickle(pickle_buffer)
        #s3.Object(s3_bucket_name, pkl_filename).put(Body=pickle_buffer.getvalue())
        
        # for local disk:
        pkl_filename = "./ceres-infer/gw_predictions/" + user_id + "_genome_wide_prediction.pkl"
        df_pred.to_pickle(pkl_filename)
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
              Input('df_pred-store-tmp-3', 'data'),
              prevent_initial_call = True) #check
def update_df_pred_store(data1, data2, pkl_filename):
    trigger_id = ctx.triggered_id
    if trigger_id=='df_pred-store-tmp-1':
        return data1
    elif trigger_id=='df_pred-store-tmp-2':
        return data2
    elif trigger_id=='df_pred-store-tmp-3':
        logging.warning('****')
        # for s3:
        #obj = s3.Object(s3_bucket_name, pkl_filename)
        #data = obj.get()['Body'].read()
        #df_pred = pickle.load(io.BytesIO(data))                 
        
        # for local:
        with open(pkl_filename , 'rb') as f:
            df_pred = pickle.load(f)
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
                                                          'font_size': '12px',
                                                          #'paddingLeft': '20px',
                                                          #'padding': '5px',
                                                          #'padding-left': '20px',
                                                            'height': '10px',
                                                            'minHeight': '10px',
                                                            'maxHeight': '10px',
                                                          },
                                              style_header={
                                                            'backgroundColor': 'dodgerblue',
                                                            'font_family': 'sans-serif',
                                                            'font_size': '12px',
                                                            #'paddingLeft': '5px',
                                                            #'padding': '5px',
                                                            'color': 'white',
                                                            'fontWeight': 'bold',
                                                            #'padding-left': '20px',
                                                            'height': '10px',
                                                            'minHeight': '10px',
                                                            'maxHeight': '10px',
                                                            'whiteSpace': 'normal',
                                                            'maxWidth': '100px',
                                                        },
                                              style_table={
                                                    'height': 400,
                                                    'overflow': 'hidden',
                                                    'textOverflow': 'ellipsis',
                                              },
                                              style_data={
                                                    'width': '50px', 'minWidth': '50px', 'maxWidth': '150px',
                                                    'overflow': 'hidden',
                                                    'textOverflow': 'ellipsis',
                                                    'minWidth': 50,
                                                    'whiteSpace': 'normal',
                                                    'height': '10px',
                                                    'minHeight': '10px',
                                                    'maxHeight': '10px',
                                                    'lineHeight': '10px'
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
              Input('distribution_color', 'value'),
              State('df_pred-store', 'data'),
              prevent_initial_call = True)
def update_overview(experiment, categories, cmap, pred_data):
    if experiment is None or pred_data is None or categories==[] or cmap is None:
        return dash.no_update
    df_pred = pd.DataFrame(pred_data)
    
    cmaph = plt.cm.get_cmap(cmap, 3)
    cmaph = [(cmaph(i)[0], cmaph(i)[1], cmaph(i)[2], cmaph(i)[3])  for i in range(3)]
    cmaph = [mpl.colors.rgb2hex(i, keep_alpha=False) for i in cmaph]

    color_dict = dict(zip(["common essential", "common nonessential", "conditional essential"], cmaph))

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
    fig.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5, tickfont_family="Arial Black",
                    showline=True, linewidth=1, linecolor='black', mirror=True,)
    fig.update_xaxes(tickfont_family="Arial Black",
                    showline=True, linewidth=1, linecolor='black', mirror=True,)
    ymax = 0

    for i, cat in enumerate(categories):
        hist_graph = go.Histogram(x=df_pred[df_pred['gene_category']==cat][experiment + ' (CERES Pred)'].values, 
                                  name=cat, 
                                  marker_color=color_dict[cat],
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
                                 line=dict(color=adjust_lightness(color_dict[cat], 0.9), width=4, dash='dot'),
                                 name=cat + ' (average)',
                                 showlegend=False))
    fig.update_layout(
    barmode='overlay',
    title_text = '<b>Distribution of Predicted CERES Scores</b>',
    title_x = 0.5,
    font_family = 'Arial',
    font_color = 'black',
    font_size=16,
    xaxis_title="<b>Predicted CERES Score</b>",
    yaxis_title="<b>Number of Genes</b>",
    legend = dict(font=dict(size=14, family="Arial",), orientation='h', yanchor='top', y=-0.2),
    margin=dict(t=50, b=80)
    )
    
    fig2 = go.Figure()
    fig2.update_layout(
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        range=[-1000, 19000],
        showticklabels=False,
        tickfont_family="Arial Black",
        showline=True, linewidth=1, linecolor='black', mirror=True,),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        tickfont_family="Arial Black",
        showline=True, linewidth=1, linecolor='black', mirror=True,))
    fig2.update_layout(plot_bgcolor='white',
    font_family = 'Arial',
    font_color = 'black',
    font_size=16,)
    fig2.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5)

    for i, cat in enumerate(categories):
        n_hist_graph = go.Scatter(x=np.sort(df_pred.sort_values(by=[experiment + ' (Z-Score)']).reindex().query('gene_category == @cat').index), 
                                  y=np.sort(df_pred[df_pred['gene_category']==cat][experiment + ' (Z-Score)'].values), 
                                  name=cat,
                                  customdata=list(df_pred[df_pred['gene_category']==cat].sort_values(by=[experiment + ' (Z-Score)']).gene),
                                  hovertemplate='gene:%{customdata}',
                                  marker_color=color_dict[cat],
                                  line=dict(color=adjust_lightness(color_dict[cat], 0.9), width=4),
                                  opacity=0.90)
        fig2.add_trace(n_hist_graph)

    fig2.update_layout(barmode='overlay')
    fig2.update_layout(
    title_text = '<b>Z-Scores of CERES Predictions</b>',
    title_x = 0.5,
    font_family = 'Arial',
    font_color = 'black',
    font_size=16,
    yaxis_title="<b>Predicted Z-Score</b>",
    #width=600,
    legend = dict(font=dict(size=14, family="Arial",), orientation='h', yanchor='top', y=-0.2),
    margin=dict(t=50, b=80)
    )

    fig3 = ff.create_distplot([df_pred[df_pred['gene_category']==cat][experiment + ' (CERES Pred)'].values for cat in categories],
                              categories,
                              colors=[color_dict[cat] for cat in categories],
                              bin_size = 0.1,
                              show_curve = True,
                              show_rug = True)
    fig3.update_layout(
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        tickfont_family="Arial Black",
        showline=True, linewidth=1, linecolor='black', mirror=True,),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        tickfont_family="Arial Black",
        showline=True, linewidth=1, linecolor='black', mirror=True,),
    xaxis1=dict(showline=True, linewidth=1, linecolor='black', mirror=True,),
    yaxis2=dict(showline=True, linewidth=1, linecolor='black', mirror=True,))
    fig3.update_layout(plot_bgcolor='white')

    fig3.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5)
    fig3.update_layout(barmode='overlay')
    fig3.update_layout(
    title_text = '<b>PDF of Predicted CERES Scores</b>',
    title_x = 0.5,
    font_family = 'Arial',
    font_color = 'black',
    font_size=16,
    xaxis_title="<b>Predicted CERES Score</b>",
    yaxis_title="<b>Probability Density</b>",  
    legend = dict(font=dict(size=14, family="Arial"), orientation='h', yanchor='top', y=-0.2),
    margin=dict(t=50, b=80)
    )
    
    fig4 = go.Figure()
    fig4.update_layout(barmode='overlay')
    fig4.update_layout(
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        tickfont_family="Arial Black",
        showline=True, linewidth=1, linecolor='black', mirror=True,),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        tickfont_family="Arial Black",
        showline=True, linewidth=1, linecolor='black', mirror=True,))
    fig4.update_layout(plot_bgcolor='white')
    fig4.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5)
    
    df_pred['P-Value'] = df_pred[experiment + ' (Z-Score)'].apply(lambda x: norm.sf(abs(x)))
    df_pred['Z-Z-Score'] = df_pred[experiment + ' (Z-Score)'].sub(df_pred[experiment + ' (Z-Score)'].mean()).div(df_pred[experiment + ' (Z-Score)'].std(ddof=0))
    df_pred['P-Value-Z'] = df_pred['Z-Z-Score'].apply(lambda x: norm.sf(abs(x)))
    df_pred['P-Value-Comb'] = df_pred['P-Value'] * df_pred['P-Value-Z']
    df_pred['log10Pcomb'] = df_pred['P-Value-Comb'].apply(lambda x: -np.log10(x))
        
    for i, cat in enumerate(categories):
        n_hist_graph_2 = go.Scatter(x=df_pred[df_pred['gene_category']==cat][experiment + ' (CERES Pred)'].values, 
                                  y=df_pred[df_pred['gene_category']==cat]['log10Pcomb'].values, 
                                  name=cat,
                                  customdata=list(df_pred[df_pred['gene_category']==cat].gene),
                                  hovertemplate='gene:%{customdata}',
                                  mode='markers',
                                  marker_color=color_dict[cat],
                                  opacity=0.90)
        fig4.add_trace(n_hist_graph_2)

    
    fig4.update_layout(
    barmode='overlay',
    title_text = '<b>Volcano Plot of CERES Scores</b>',
    title_x = 0.5,
    font_family = 'Arial',
    font_color = 'black',
    font_size=16,
    xaxis_title="<b>Predicted CERES Score</b>",
    yaxis_title="<b>-log10(p-value)</b>",
    legend = dict(font=dict(size=14, family="Arial"), orientation='h', yanchor='top', y=-0.2),
    margin=dict(t=50, b=80)
    )

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

            return dbc.Table.from_dataframe(df_pred, striped=True, bordered=True, hover=True, size='sm')
        
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
            return dbc.Table.from_dataframe(df_pred, striped=True, bordered=True, hover=True, size='sm')

@app.callback(Output('output_z_score', 'figure'),
              Output('selective-gene-table', 'children'),
              Output('download-z-genes-button', 'disabled'),
              Output('run-gprofiler-zscores', 'disabled'),
              Output('z-score-hits', 'data'),
              Input('choose_experiment_3', 'value'),
              State('experiment_labels', 'data'),
              State('df_pred-store', 'data'),
              Input('select_number_genes', 'value'),
              Input('primary_c', 'value'),
              Input('experiment_c', 'value'),
              Input('sig_genes_threshold_1', 'value'),
              Input('z-score_threshold_1', 'value'),
              State('select_primary_cat', 'value'),
              Input('select_sub_cat', 'value'),
              Input('secondary_c', 'value'),
              prevent_initial_call = True)
def update_z_score_graph(experiment, experiments, pred_data, n_genes, primary_c, experiment_c, sig_thres, avg_range, cat, sub_cat, secondary_c):
    if experiment is None or experiments is None or pred_data is None or sig_thres is None or avg_range is None:
        return dash.no_update
    else:
        df_pred = pd.DataFrame(pred_data)
        df_pred = df_pred[(df_pred['avg'] >= avg_range[0]) & (df_pred['avg'] <= avg_range[1]) & (df_pred[experiment + ' (CERES Pred)'] <= sig_thres)]
        
        color_list = ['#7f00ff', '#0000ff', '#00ffff','#ffff00', '#007fff', '#7fff00', '#00ff00']
        color_dict = dict(zip(experiments, color_list))

        fig = go.Figure()

        df_pred = df_pred.sort_values(by=experiment + ' (Z-Score)')
        top_genes = df_pred['gene'][0:int(n_genes)]
        n_genes = len(top_genes)

        plot_min = 0
        plot_max = 0

        dm_data_small = dm_data.dropna()
        dm_data_small = dm_data_small.merge(df_sample_data, how='left', right_on='DepMap_ID', left_index=True)
        
        fig.add_shape(type="line",
            x0=-1,
            x1=-1,
            y0=-1,
            y1=n_genes-0.7,
            line=dict(color='red', 
                        width=3,
                        dash='dot'))
        
        fig.add_annotation(x=-1, y=n_genes,
        text="<b>Avg CERES Score:</b><br>Common Essential",
        showarrow=False,
        #xshift=,
        #yshift=-1,
        font=dict(size=12,
                  color='red'))
        
        fig.add_shape(type="line",
            x0=0,
            x1=0,
            y0=-1,
            y1=n_genes-0.7,
            line=dict(color='blue', 
                        width=3,
                        dash='dot'))
        
        fig.add_annotation(x=0, y=n_genes,
            text="<b>Avg CERES Score:</b><br>Common Nonessential",
            showarrow=False,
            #xshift=,
            #yshift=-1,
            font=dict(size=12,
                  color='blue'))

        for i, gene in reversed(list(enumerate(top_genes))):
            if i > 0:
                show = False
            else:
                show = True
            
            x = dm_data_small[gene].values[~np.isnan(dm_data_small[gene].values)]

            if min(x) < plot_min:
                plot_min = min(x)
            if max(x) > plot_max:
                plot_max = max(x)

            fig.add_trace(go.Violin(x=x,
                                    y=[gene] * len(x),
                                    side='positive',
                                    line_color=primary_c,
                                    name='DepMap: All Cells',
                                    orientation='h',
                                    showlegend=show,
                                    points=False,
                                    meanline_visible=True))
            
            if cat is not None and sub_cat is not None:
                y = dm_data_small[dm_data_small[cat]==sub_cat][gene].values[~np.isnan(dm_data_small[dm_data_small[cat]==sub_cat][gene])]

                if min(y) < plot_min:
                    plot_min = min(y)
                if max(y) > plot_max:
                    plot_max = max(y)

                fig.add_trace(go.Violin(x=y,
                                y=[gene] * len(y),
                                side='negative',
                                line_color=secondary_c,
                                name='DepMap: ' + sub_cat,
                                orientation='h',
                                showlegend=show,
                                points=False,
                                meanline_visible=True)) 

            fig.add_trace(go.Scatter(x=df_pred[experiment + ' (CERES Pred)'][df_pred['gene']==gene],
                                     y=[gene],
                                     mode='markers',
                                     marker=dict(
                                         symbol='diamond',
                                         color=experiment_c,
                                         size=9,
                                         line=dict(
                                             color='DarkSlateGray',
                                             width=1
                                         )), name=experiment + ' (CERES Prediction)',
                                         showlegend=show))

        fig.add_trace(go.Scatter(x=None, y=None, name="ghost1", showlegend=False))
        fig.update_layout(xaxis2= {'anchor': 'y', 'overlaying': 'x', 'side': 'top'},
                  yaxis_domain=[0, 1])

        if n_genes > 20:
            height = 30 * n_genes
        else:
            height = 610
            
        leg_y = -0.1 * 20 / n_genes
        
        fig.update_layout(
        plot_bgcolor='white',
        title_text = '<b>Top {} Selective Essential Genes by Z-Score</b><br><sup>{} (CERES Prediction)</sup>'.format(n_genes, experiment),
        title_x = 0.5,
        #title_y = 1,
        xaxis_title="<b>CERES Score</b>",
        font_family = 'Arial',
        font_color = 'black',
        font_size=16,
        legend = dict(font=dict(size=16, family="Arial Black", color="black"), orientation='h', yanchor='top', x=-0.1, y=leg_y),
        margin=dict(t=80, b=0),
        violingap=0,
        violingroupgap=0,
        violinmode='overlay',
        height=height,
        #width=600,
        )

        fig.update_yaxes(
            range=[-1, n_genes],
            zerolinecolor='black',
            zerolinewidth=0.5,
            tickfont_family="Arial Black", 
            tickfont=dict(size=16)
        )
        

        fig.update_xaxes(dict(
            zeroline=True,
            zerolinecolor='lightgray',
            zerolinewidth=0.5,
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray'),
            range=[plot_min-0.2, plot_max+0.2],
            tickfont_family="Arial Black", 
            tickfont=dict(size=16))
        
        minx = round(fig.layout.xaxis.range[0], 1)
        maxx = round(fig.layout.xaxis.range[1], 1)

        x_tick_vals = [i/10 for i in range(int(10*minx), int(10*maxx)) if i%5==0]
        
        s_c = {-4.0: ['black', '16px'],
               -3.5: ['black', '16px'],
               -3.0: ['black', '16px'],
               -2.5: ['black', '16px'],
               -2.0: ['black', '16px'],
               -1.5: ['black', '16px'],
               -1.0: ['red', '18px'],
               -0.5: ['black', '16px'],
                0.0: ['blue', '18px'],
                0.5: ['black', '16px'],
                1.0: ['black', '16px'],
                1.5: ['black', '16px'],
                2.0: ['black', '16px'],
                2.5: ['black', '16px']}
        
        colors = [s_c[i] for i in x_tick_vals if i in s_c.keys()]
        keys = dict(zip(x_tick_vals, colors))
        ticktext = [color(v, k) for k, v in keys.items()]
        
        fig.update_layout(
            xaxis=dict(tickmode='array', ticktext=ticktext, tickvals=x_tick_vals),
        )
        
        fig.add_shape(type="line",
            x0=plot_min-0.2,
            x1=plot_max+0.2,
            y0=-1,
            y1=-1,
            line=dict(color='black', 
                        width=4,))

        df_return = df_pred[['gene', experiment + ' (CERES Pred)', experiment + ' (Z-Score)']][0:int(n_genes)]
        df_return = df_return.rename(columns={'gene': 'Gene', experiment + ' (CERES Pred)': 'CERES Pred', experiment + ' (Z-Score)': 'Z-Score'})
        df_return = df_return.round(3)
        select_genes = df_return['Gene'].values.tolist()
        df_return['Gene'] = df_return['Gene'].apply(lambda x: create_gene_link(x, 'depmap', 'gene'))

        return fig, dbc.Table.from_dataframe(df_return, striped=True, bordered=True, hover=True, size='sm'), False, False, select_genes

@app.callback(Output('download-z-genes', 'data'),
              State('z-score-hits', 'data'),
              State('choose_experiment_3', 'value'),
              Input('download-z-genes-button', 'n_clicks'),
              prevent_initial_call = True) #check
def download_z_genes(genes, experiment, n_clicks):
    if genes is None or experiment is None:
        return dash.no_update
    else:
        return_genes = "\n".join(genes)
        return dict(content=return_genes, filename=experiment+'_selective_genes.txt')

@app.callback(Output('download-p-genes', 'data'),
              State('p-value-hits', 'data'),
              State('choose_experiment_4', 'value'),
              Input('download-p-genes-button', 'n_clicks'),
              prevent_initial_call = True) #check
def download_z_genes(genes, experiment, n_clicks):
    if genes is None or experiment is None:
        return dash.no_update
    else:
        return_genes = "\n".join(genes)
        return dict(content=return_genes, filename=experiment+'_selective_genes.txt')

@app.callback(Output('enrichment-hits', 'data'),
              Input('p-value-hits', 'data'),
              Input('z-score-hits', 'data'),
              prevent_initial_call = True)
def update_enrichment(p_hits, z_hits):
    if p_hits is None and z_hits is None:
        return dash.no_update
    else:
        trigger_id = ctx.triggered_id
        if trigger_id == 'p-value-hits':
            return p_hits
        elif trigger_id == 'z-score-hits':
            return z_hits
        

@app.callback(Output('url', 'href'),
              Input('continue-gprofile', 'n_clicks'),
              prevent_initial_call = True) #check
def redirect_to_gprofiler_1(n_clicks):
    if n_clicks:
        return "/enrichment"

@app.callback(Output('manhattan-plot', 'figure'),
              Output('output-term-table', 'children'),
              Output('gprofile-gene-mapping', 'children'),
              Output('meta-data-info', 'children'),
              State('enrichment-hits', 'data'),
              Input('dummy_div_enrichment', 'children'),)
def generate_manhattan_plot(genes, aux):
    if genes is None:
        return dash.no_update
    else:
        r = requests.post(
            url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
            json={
                'organism':'hsapiens',
                'query': genes,
                'highlight': True,
            })
        fig = gostplot(r.json())
        return fig   

@app.callback(Output('select_sub_cat', 'options'),
              Input('select_primary_cat', 'value'),
              prevent_initial_call = True)
def return_options(color_by):
    dm_data_small = dm_data.dropna()
    dm_data_small = dm_data_small.merge(df_sample_data, how='left', right_on='DepMap_ID', left_index=True)
    values = dm_data_small[color_by].unique()
    return [{'label': value, 'value': value} for value in values]

@app.callback(Output("m_select_sub_cat", 'options'),
              Output("m_select_sub_cat", 'value'),
              Input("m_select_secondary_cat", 'value'),
              prevent_initial_call = True)
def return_options(color_by):
    dm_data_small = dm_data.dropna()
    dm_data_small = dm_data_small.merge(df_sample_data, how='left', right_on='DepMap_ID', left_index=True)
    values = dm_data_small[color_by].unique()
    values = [i for i in values if i == i]
    return [{'label': value, 'value': value} for value in values], None


@app.callback(Output('output_volcano_plot', 'figure'),
              Output('significant-gene-table-p', 'children'),
              Output('download-p-genes-button', 'disabled'),
              Output('run-gprofiler-pvalue', 'disabled'),
              Output('p-value-hits', 'data'),
              Input('choose_experiment_4', 'value'),
              Input('choose_gene_cats', 'value'),
              Input('conditions_c', 'value'),
              Input('significant_c', 'value'),
              Input('gene_effect_c', 'value'),
              Input('p_val_c', 'value'),
              Input('ceres_range_slider', 'value'),
              Input('p_val_threshold', 'value'),
              State('df_pred-store', 'data'),
              prevent_initial_call = True)
def update_volcano(experiment, categories, cmap, c_sig, c_ge, c_p, ceres_thres, p_val_thres, pred_data):
    if experiment is None or pred_data is None or categories==[] or categories is None or cmap is None or c_sig is None or c_ge is None or c_p is None or ceres_thres==[] or p_val_thres is None:
        return dash.no_update
    df_pred = pd.DataFrame(pred_data)
    
    cmaph = plt.cm.get_cmap(cmap, 3)
    cmaph = [(cmaph(i)[0], cmaph(i)[1], cmaph(i)[2], cmaph(i)[3])  for i in range(3)]
    cmaph = [mpl.colors.rgb2hex(i, keep_alpha=False) for i in cmaph]

    color_dict = dict(zip(["common essential", "common nonessential", "conditional essential"], cmaph))
    
    df_pred['P-Value'] = df_pred[experiment + ' (Z-Score)'].apply(lambda x: norm.sf(abs(x)))
    df_pred['Z-Z-Score'] = df_pred[experiment + ' (Z-Score)'].sub(df_pred[experiment + ' (Z-Score)'].mean()).div(df_pred[experiment + ' (Z-Score)'].std(ddof=0))
    df_pred['P-Value-Z'] = df_pred['Z-Z-Score'].apply(lambda x: norm.sf(abs(x)))
    df_pred['P-Value-Comb'] = df_pred['P-Value'] * df_pred['P-Value-Z']
    df_pred['log10Pcomb'] = df_pred['P-Value-Comb'].apply(lambda x: -np.log10(x))
    
    fig = go.Figure()
    fig.update_layout(barmode='overlay')
    fig.update_layout(
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        tickfont_family="Arial Black",
        showline=True, linewidth=1, linecolor='black', mirror=True,),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        tickfont_family="Arial Black",
        showline=True, linewidth=1, linecolor='black', mirror=True,))
    fig.update_layout(plot_bgcolor='white')
    fig.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=1.5)
        
    for i, cat in enumerate(categories):
        n_hist_graph_2 = go.Scatter(x=df_pred[(df_pred['gene_category']==cat) & (((df_pred[experiment + ' (CERES Pred)'] > ceres_thres[0]) & (df_pred[experiment + ' (CERES Pred)'] < ceres_thres[1])) | (df_pred['log10Pcomb'] < p_val_thres))][experiment + ' (CERES Pred)'].values, 
                                  y=df_pred[(df_pred['gene_category']==cat) & (((df_pred[experiment + ' (CERES Pred)'] > ceres_thres[0]) & (df_pred[experiment + ' (CERES Pred)'] < ceres_thres[1])) | (df_pred['log10Pcomb'] < p_val_thres))]['log10Pcomb'].values, 
                                  name=cat,
                                  customdata=list(df_pred[(df_pred['gene_category']==cat) & (((df_pred[experiment + ' (CERES Pred)'] > ceres_thres[0]) & (df_pred[experiment + ' (CERES Pred)'] < ceres_thres[1])) | (df_pred['log10Pcomb'] < p_val_thres))].gene),
                                  hovertemplate='gene:%{customdata}',
                                  mode='markers',
                                  marker_color=color_dict[cat],
                                  #line=dict(color=adjust_lightness(all_colors[color_dict[cat]], 0.9), width=4),
                                  opacity=0.30)
        fig.add_trace(n_hist_graph_2)
        
    fig.add_trace(
        go.Scatter(x=df_pred[(df_pred['gene_category'].isin(categories)) & ~(((df_pred[experiment + ' (CERES Pred)'] > ceres_thres[0]) & (df_pred[experiment + ' (CERES Pred)'] < ceres_thres[1])) | (df_pred['log10Pcomb'] < p_val_thres))][experiment + ' (CERES Pred)'].values, 
                                  y=df_pred[(df_pred['gene_category'].isin(categories)) & ~(((df_pred[experiment + ' (CERES Pred)'] > ceres_thres[0]) & (df_pred[experiment + ' (CERES Pred)'] < ceres_thres[1])) | (df_pred['log10Pcomb'] < p_val_thres))]['log10Pcomb'].values, 
                                  name='significant genes',
                                  customdata=list(df_pred[(df_pred['gene_category'].isin(categories)) & ~(((df_pred[experiment + ' (CERES Pred)'] > ceres_thres[0]) & (df_pred[experiment + ' (CERES Pred)'] < ceres_thres[1])) | (df_pred['log10Pcomb'] < p_val_thres))].gene),
                                  hovertemplate='gene:%{customdata}',
                                  mode='markers',
                                  marker_color=c_sig,
                                  #line=dict(color=adjust_lightness(all_colors[color_dict[cat]], 0.9), width=4),
                                  opacity=0.80)
    )
    
    fig.add_trace(
        go.Scatter(x=[ceres_thres[0], ceres_thres[0]], 
                   y=[min(df_pred['log10Pcomb'].values), max(df_pred['log10Pcomb'].values)], 
                   name='min effect size',
                   mode='lines',
                   line = dict(color=c_ge, width=3, dash='dash'),
                  showlegend =False)
    )
    
    fig.add_trace(
        go.Scatter(x=[ceres_thres[1], ceres_thres[1]], 
                   y=[min(df_pred['log10Pcomb'].values), max(df_pred['log10Pcomb'].values)], 
                   name='max effect size',
                   mode='lines',
                   line = dict(color=c_ge, width=3, dash='dash'),
                  showlegend=False)
    )

    fig.add_trace(
        go.Scatter(x=[min(df_pred[experiment + ' (CERES Pred)']), max(df_pred[experiment + ' (CERES Pred)'])], 
                   y=[p_val_thres, p_val_thres],
                   name='min p-value',
                   mode='lines',
                   line = dict(color=c_p, width=3, dash='dash'),
                  showlegend=False)
    )

    
    fig.update_layout(
    barmode='overlay',
    title_text = '<b>Volcano Plot of CERES Scores</b>',
    title_x = 0.5,
    font_family = 'Arial',
    font_color = 'black',
    font_size=16,
    xaxis_title="<b>Predicted CERES Score</b>",
    yaxis_title="<b>-log10(p-value)</b>",
    #width=600,
    height=700,
    legend = dict(font=dict(size=16, family="Arial"), orientation='h', yanchor='top', y=-0.2),
    #margin=dict(t=50, b=80)
    )

    df_return = df_pred[(df_pred['gene_category'].isin(categories)) & ~(((df_pred[experiment + ' (CERES Pred)'] > ceres_thres[0]) & (df_pred[experiment + ' (CERES Pred)'] < ceres_thres[1])) | (df_pred['log10Pcomb'] < p_val_thres))][['gene', experiment + ' (CERES Pred)', experiment + ' (Z-Score)']]
    df_return = df_return.rename(columns={'gene': 'Gene', experiment + ' (CERES Pred)': 'CERES Pred', experiment + ' (Z-Score)': 'Z-Score'})
    df_return = df_return.round(3)
    select_genes = df_return['Gene'].values.tolist()
    df_return['Gene'] = df_return['Gene'].apply(lambda x: create_gene_link(x, 'depmap', 'gene'))

    return fig, dbc.Table.from_dataframe(df_return, striped=True, bordered=True, hover=True, size='sm'), False, False, select_genes

#-----------------------------------------------------------GENES PAGE CALLBACKS------------------------------------------------------------------------------

@app.callback(Output('pick-community-1', 'options'),
          Output('pick-community-2', 'options'),
          Output('pick-community-3', 'options'),
          Input('dummy_div_genes', 'children'),)
def update_dropdown(aux):

    landmark_genes = [{'label': 'louvain_' + str(i), 'value': i} for i in range(44)]
    landmark_genes.insert(0, {'label': 'ALL', 'value': 'ALL'})

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
            clust = modularity[community]
            
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
            clust = modularity[community]
            
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
            clust = modularity[community]
            
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
    fig.update_layout(
    title_text = '<b>Distribution of all CERES Scores for ' + gene2visualize + '</b>',
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
    fig.update_layout(
    title_text = '<b>Distribution of all CERES Scores for ' + gene2visualize + '</b>',
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
            title_text = '<b>Pairwise CERES Score comparison between ' + gene_1 + ' and ' + gene_2 + '</b>',
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

        fig.update_layout(height=550)

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
                                         size=5,
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
            title_text = '<b>Pairwise CERES Score comparison between ' + gene_1 + ' and ' + gene_2 + '</b>',
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

        fig.update_layout(height=550)

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
              Output('multi_gene_comp_network', 'children'),
              Input('multi-gene-list-dropdown', 'value'),
              Input('m_primary_color', 'value'),
              Input('m_experiment_color', 'value'),
              State('m_select_secondary_cat', 'value'),
              Input('m_select_sub_cat', 'value'),
              Input('m_secondary_color', 'value'),
              State('experiment_labels', 'data'),
              State('df_pred-store', 'data'),
              prevent_inital_call=True)
def run_multiple_genes(gene_list, primary_c, cmap, cat, sub_cat, secondary_c, experiments, pred_data):
    if gene_list is None or experiments is None or pred_data is None or primary_c is None or cmap is None or cat is None or sub_cat is None or secondary_c is None:
        raise PreventUpdate
    elif len(gene_list) < 2 or len(gene_list) > 6:
        raise PreventUpdate
    else:
        n = len(gene_list)
        fig = make_subplots(rows=n, cols=n, shared_yaxes=True, shared_xaxes=True, horizontal_spacing=0.01, vertical_spacing=0.01)
        df_pred = pd.DataFrame(pred_data)
        dm_data_small = dm_data.dropna()
        dm_data_small = dm_data_small.merge(df_sample_data, how='left', right_on='DepMap_ID', left_index=True)
    
        cmaph = plt.cm.get_cmap(cmap, len(experiments))
        cmaph = [(cmaph(i)[0], cmaph(i)[1], cmaph(i)[2], cmaph(i)[3])  for i in range(len(experiments))]
        cmaph = [mpl.colors.rgb2hex(i, keep_alpha=False) for i in cmaph]
        color_dict = dict(zip(experiments, cmaph))
        
        symbol_list = ['diamond', 'diamond', 'diamond']
        symbol_dict = dict(zip(experiments, iter(symbol_list)))
    
        for i in range (1,n+1):
            fig.update_yaxes(title_text="<b>"+gene_list[i-1]+"</b>", row=i, col=1)
            for j in range (1,n+1):
                if gene_list[j-1] == gene_list[i-1]:
                    a = dm_data_small[gene_list[j-1]]
                    
                    Q1a = np.percentile(a, 25, method='midpoint')
                    Q3a = np.percentile(a, 75, method='midpoint')
                    IQRa = Q3a - Q1a
                    
                    uppera=Q3a+1.5*IQRa
                    lowera=Q1a-1.5*IQRa
                    
                    fig2 = ff.create_distplot([dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                             (dm_data_small[gene_list[j-1]]<=uppera)][gene_list[j-1]].values],
                                              [gene_list[j-1]],
                                              histnorm='probability density',
                                              bin_size = 0.1,
                                              show_curve = True,
                                              show_rug = False)
                    
                    fig.add_trace(go.Histogram(fig2['data'][0],
                           marker_color='cyan', showlegend=False), row=i, col=j)
                    
                    fig.add_trace(go.Scatter(fig2['data'][1],
                         line=dict(color='blue', width=4),
                        showlegend=False), row=i, col=j)
                    
                    if i == n:
                        fig.update_xaxes(title_text="<b>"+gene_list[j-1]+"</b>", row=i, col=j)
                    
                if not gene_list[j-1] == gene_list[i-1]:
                    
                    rg=LinearRegression()
                    a = dm_data_small[gene_list[j-1]]
                    b = dm_data_small[gene_list[i-1]]

                    Q1a = np.percentile(a, 25, method='midpoint')
                    Q3a = np.percentile(a, 75, method='midpoint')
                    IQRa = Q3a - Q1a

                    Q1b = np.percentile(b, 25, method='midpoint')
                    Q3b = np.percentile(b, 75, method='midpoint')
                    IQRb = Q3b - Q1b

                    uppera=Q3a+1.5*IQRa
                    lowera=Q1a-1.5*IQRa

                    upperb=Q3b+1.5*IQRb
                    lowerb=Q1b-1.5*IQRb
                    
                    if i > j:
                        if i == 2:
                            show = True
                        else:
                            show = False
                            
                        r, _ = pearsonr(dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                                 (dm_data_small[gene_list[j-1]]<=uppera) &
                                                                 (dm_data_small[gene_list[i-1]]>=lowerb) &
                                                                 (dm_data_small[gene_list[i-1]]<=upperb)][gene_list[j-1]], 
                                        dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                                 (dm_data_small[gene_list[j-1]]<=uppera) &
                                                                 (dm_data_small[gene_list[i-1]]>=lowerb) &
                                                                 (dm_data_small[gene_list[i-1]]<=upperb)][gene_list[i-1]])

                        fig.add_trace(go.Scatter(x=dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                                 (dm_data_small[gene_list[j-1]]<=uppera) &
                                                                 (dm_data_small[gene_list[i-1]]>=lowerb) &
                                                                 (dm_data_small[gene_list[i-1]]<=upperb)
                                                                ][gene_list[j-1]], 
                                                 y=dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                                 (dm_data_small[gene_list[j-1]]<=uppera) &
                                                                 (dm_data_small[gene_list[i-1]]>=lowerb) &
                                                                 (dm_data_small[gene_list[i-1]]<=upperb)][gene_list[i-1]],
                                                 text=dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                                 (dm_data_small[gene_list[j-1]]<=uppera) &
                                                                 (dm_data_small[gene_list[i-1]]>=lowerb) &
                                                                 (dm_data_small[gene_list[i-1]]<=upperb)]['stripped_cell_line_name'],
                                                 mode='markers',
                        marker=dict(
                        size=10,
                        color = primary_c,
                        opacity = 0.4,line=dict(
                            color=primary_c,
                            width=0.2
                        )),
                        name = "DepMap: All Cells",
                        showlegend=show,
                        ), 
                        row=i, col=j)
                        
                        fig.add_annotation(text="<b>r = {0:.2f}</b>".format(r),
                              xref="paper", yref="paper",
                              x=Q3a, 
                              y=Q3b+2.5*IQRb, 
                              showarrow=False, 
                              font=dict(size=23, color='red'), 
                              row=i, col=j)

                        for e in experiments:

                            gene_1_pred = df_pred[df_pred['gene'] == [gene_list[j-1]][0]][e + ' (CERES Pred)'].values
                            gene_2_pred = df_pred[df_pred['gene'] == [gene_list[i-1]][0]][e + ' (CERES Pred)'].values

                            fig.add_trace(go.Scatter(x=gene_1_pred, y=gene_2_pred, mode='markers', 
                            marker=dict(
                            symbol=symbol_dict[e],
                            color=color_dict[e],
                            size=10,
                            opacity = 1,
                            line=dict(
                                color='black',
                                width=1
                            )), name=e + ' (CERES Prediction)',
                            showlegend=show), row=i, col=j)
                            
                        if i == n:
                            fig.update_xaxes(title_text="<b>"+gene_list[j-1]+"</b>", row=i, col=j)
                            
                    if i < j:
                        if j == 2:
                            show = True
                        else:
                            show = False

                        r, _ = pearsonr(dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                                 (dm_data_small[gene_list[j-1]]<=uppera) &
                                                                 (dm_data_small[gene_list[i-1]]>=lowerb) &
                                                                 (dm_data_small[gene_list[i-1]]<=upperb) &
                                                                 (dm_data_small[cat]==sub_cat)][gene_list[j-1]], 
                                        dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                                 (dm_data_small[gene_list[j-1]]<=uppera) &
                                                                 (dm_data_small[gene_list[i-1]]>=lowerb) &
                                                                 (dm_data_small[gene_list[i-1]]<=upperb) &
                                                                 (dm_data_small[cat]==sub_cat)][gene_list[i-1]])
                        
                        fig.add_trace(go.Scatter(x=dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                                 (dm_data_small[gene_list[j-1]]<=uppera) &
                                                                 (dm_data_small[gene_list[i-1]]>=lowerb) &
                                                                 (dm_data_small[gene_list[i-1]]<=upperb) &
                                                                 (dm_data_small[cat]==sub_cat)][gene_list[j-1]], 
                                                 y=dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                                 (dm_data_small[gene_list[j-1]]<=uppera) &
                                                                 (dm_data_small[gene_list[i-1]]>=lowerb) &
                                                                 (dm_data_small[gene_list[i-1]]<=upperb) &
                                                                 (dm_data_small[cat]==sub_cat)][gene_list[i-1]],
                                                 text=dm_data_small[(dm_data_small[gene_list[j-1]]>=lowera) &
                                                                 (dm_data_small[gene_list[j-1]]<=uppera) &
                                                                 (dm_data_small[gene_list[i-1]]>=lowerb) &
                                                                 (dm_data_small[gene_list[i-1]]<=upperb) &
                                                                 (dm_data_small[cat]==sub_cat)]['stripped_cell_line_name'],
                                                 mode='markers',
                        marker=dict(
                        size=10,
                        color = secondary_c,
                        opacity = 0.4,line=dict(
                            color=secondary_c,
                            width=0.2
                        )),
                        name = "DepMap: " + sub_cat,
                        showlegend=show,
                        ), row=i, col=j)

                        fig.add_annotation(text="<b>r = {0:.2f}</b>".format(r),
                            xref="paper", yref="paper",
                            x=Q3a, 
                            y=Q3b+2.5*IQRb, 
                            showarrow=False, 
                            font=dict(size=23, color='red'), 
                            row=i, col=j)

                        for e in experiments:

                            gene_1_pred = df_pred[df_pred['gene'] == [gene_list[j-1]][0]][e + ' (CERES Pred)'].values
                            gene_2_pred = df_pred[df_pred['gene'] == [gene_list[i-1]][0]][e + ' (CERES Pred)'].values

                            fig.add_trace(go.Scatter(x=gene_1_pred, y=gene_2_pred, mode='markers', 
                            marker=dict(
                            symbol=symbol_dict[e],
                            color=color_dict[e],
                            size=10,
                            opacity = 1,
                            line=dict(
                                color='black',
                                width=1
                            )), name=e + " (CERES Prediction)",
                            showlegend = False), 
                            row=i, col=j)

                        if i == n:
                            fig.update_xaxes(title_text="<b>"+gene_list[j-1]+"</b>", row=i, col=j)

        fig.update_layout(showlegend=True,
                          plot_bgcolor='white',
                          title_text='<b>Multiple Pairwise Genes Regression Analysis</b>',
                          font_family = 'Arial',
                          font_color = 'black',
                          font_size=22,
                          title_x=0.5,
                          #legend_title="Legend",
                          legend = dict(font=dict(size=20, family="Arial", color='black'), orientation='h', yanchor='top', itemwidth=30))
        fig.update_yaxes(zeroline=True, zerolinecolor='lightgray', zerolinewidth=0.8,
                        showline=True, linewidth=1.5, linecolor='black', mirror=True,
                        tickfont_family="Arial Black", tickfont=dict(size=16))
        fig.update_xaxes(zeroline=True, zerolinecolor='lightgray', zerolinewidth=0.8,
                        showline=True, linewidth=1.5, linecolor='black', mirror=True,
                        tickfont_family="Arial Black", tickfont=dict(size=16))
        fig.update_xaxes(dict(
            showgrid=True,
            gridwidth=0.8,
            gridcolor='lightgray'),)
        fig.update_yaxes(dict(
            showgrid=True,
            gridwidth=0.8,
            gridcolor='lightgray'),)
        #fig['layout'].update(height=1000)
        
        
        clust = gene_list.copy()
        for c,j in G.edges(gene_list):
            clust.append(j)
                
        clust = list(set(clust))

        subgraph = nx.Graph()
        subgraph.add_nodes_from(clust)
        subgraph.add_edges_from(G.edges(gene_list))
        sizes_dict = nx.degree_centrality(subgraph)
        low_s = min([v for k, v in sizes_dict.items()])
        high_s = max([v for k, v in sizes_dict.items()])
        map_data = 'mapData(size, ' + str(low_s) + ', ' + str(high_s) + ', 20, 100)'

        n_colors = ['red' if i in gene_list else 'blue' for i in clust]
        n_class_dict = dict(zip(clust, n_colors))
        
        nodes = [
                {
                    'data': {'id': label, 'label': label, 'size': sizes_dict[label]},
                    'classes': n_class_dict[label],
                }
                for label in clust
            ]

        edges = [
            {
                'data': {'source': source, 'target': target},
            }
            for source, target in G.edges(gene_list)
        ]

        elements = nodes + edges

        group_select_dict = {
                        'selector': 'node',
                        'style': {
                            'content': 'data(label)',
                            'width': map_data,
                            'height': map_data,
                        }
                    }
        
        network =  [html.Br(),
                    dbc.Row([
                        dbc.Col(width=3),
                        dbc.Col(dbc.Label("Gene Network", style={'text-align': 'center', 'font-size': '200%', "font-weight": "bold", 'font-family': 'Arial'}), width=6),
                        dbc.Col(width=3),
                    ]),
                    html.Br(),
                    cyto.Cytoscape(
                    id='cytoscape-layout-2',
                    elements=elements,
                    style={'width': '90%', 'height': '700px'},
                    layout={
                        'name': 'cose',
                        'animate': True
                    },
                    stylesheet=[
                        # Group selectors
                        group_select_dict,

                        # Class selectors
                        {
                            'selector': '.red',
                            'style': {
                                'background-color': 'red',
                                'line-color': 'magenta'
                            }
                        },
                        {
                            'selector': '.blue',
                            'style': {
                                'background-color': 'lightblue',
                                'line-color': 'grey'
                            }
                        }
                    ]
                )]

        return fig, network

#----------------------------------------------------------------------Community Page Callbacks---------------------------------------


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

@app.callback(Output('network_graph_1', 'children'),
              Output('louvain-gene-table', 'children'),
              Input('select_community', 'value'),
              Input('cell-line-or-experiment', 'data'),
              State('select_experiment', 'value'),
              State('select_experiment', 'options'),
              State('select_experiment_cluster', 'value'),
              State('df_pred-store', 'data'),
              Input('select_network_layout', 'value'),
              Input('sig_genes_threshold', 'value'),
              Input('z-score_threshold', 'value'),
              Input('select_node_color_by', 'value'),
              Input('select_node_colormap', 'value'),
              Input('cres_threshold', 'value'),
              Input('select_node_size_by', 'value'),
              )
def create_network_graph_1(feature, choice, cell_line, cell_line_options, 
                           experiment, pred_data, layout, c_threshold, z_threshold, 
                           color_by, cmap, cres, size_by):
    if feature is None or (cell_line is None and experiment is None) or layout is None or color_by is None or cmap is None or size_by is None:
        return dash.no_update
    else:
        
        if color_by == 'community':
            n_colors_dict = {}        
            cmaph = plt.cm.get_cmap(cmap, len(feature))
            cmaph = [(255*cmaph(i)[0], 255*cmaph(i)[1], 255*cmaph(i)[2])  for i in range(len(feature))]
        
        clust = []
        
        for i in range(len(feature)):
            clust1 = []
            clust1.extend(modularity[i])
            clust.extend(modularity[i])
            
            if color_by == 'community':
                
                for c,j in G.edges(clust1):
                    clust1.append(j)
                    
                n_colors_dict.update({a: cmaph[i] for a in clust1})
        
        clust2 = clust.copy()
        for c,j in G.edges(clust):
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
            cell_line_label = [x['label'] for x in cell_line_options if x['value'] == cell_line][0]
            subtitle = cell_line_label

        elif choice == 'experiment':

            df_pred = pd.DataFrame(pred_data)
            df_pred = df_pred[['gene', 'gene_category', 'avg', 'std', experiment + ' (CERES Pred)', experiment + ' (Z-Score)']]
            df_pred = df_pred[df_pred['gene'].isin(clust2)]

            df_return = df_pred[['gene', experiment + ' (CERES Pred)', experiment + ' (Z-Score)']]
            
            df_return = df_return.rename(columns={'gene': 'Gene', experiment + ' (CERES Pred)': 'CERES Pred', experiment + ' (Z-Score)': 'Z-Score'})
            df_return = df_return.round(3)

            sig_genes = df_return[(df_return['CERES Pred'] < c_threshold) & (abs(df_return['Z-Score']) >= z_threshold)]['Gene'].values.tolist()
            subtitle = experiment + ' (CERES Prediction)'
        
        n_colors = ['red' if i in sig_genes else 'blue' for i in clust2]
        n_class_dict = dict(zip(clust2, n_colors))
        
        e_colors = ['red' if source in sig_genes and target in sig_genes else None for source, target in G.edges(clust)]
        e_class_dict = dict(zip(G.edges(clust), e_colors))
        
        if size_by == 'ceres_score':
            if choice == 'cell line':
                sizes = [df_return[df_return['Gene'] == i]['CERES Score'].values.tolist()[0] if i in df_return['Gene'].values.tolist() else 0 for i in clust2]
                sizes_dict = dict(zip(clust2, sizes))
            elif choice == 'experiment':
                sizes = [df_return[df_return['Gene'] == i]['CERES Pred'].values.tolist()[0] if i in df_return['Gene'].values.tolist() else 0 for i in clust2]
                sizes_dict = dict(zip(clust2, sizes))
        elif size_by == 'z_score':
            sizes = [abs(df_return[df_return['Gene'] == i]['Z-Score'].values.tolist()[0]) if i in df_return['Gene'].values.tolist() else 0 for i in clust2]
            sizes_dict = dict(zip(clust2, sizes))
        elif size_by == 'degree_centrality':
            sizes_dict = nx.degree_centrality(subgraph)
        elif size_by == 'closeness_centrality':
            sizes_dict = nx.closeness_centrality(subgraph)
        elif size_by == 'betweenness_centrality':
            sizes_dict = nx.betweenness_centrality(subgraph)
                    
        if color_by == 'ceres_score':
            if choice == 'cell line':
                n_colors = [df_return[df_return['Gene'] == i]['CERES Score'].values.tolist()[0] if i in df_return['Gene'].values.tolist() else 0 for i in clust2]
                n_colors_dict = dict(zip(clust2, n_colors))
            elif choice == 'experiment':
                n_colors = [df_return[df_return['Gene'] == i]['CERES Pred'].values.tolist()[0] if i in df_return['Gene'].values.tolist() else 0 for i in clust2]
                n_colors_dict = dict(zip(clust2, n_colors))
        elif color_by == 'z_score':
            n_colors = [abs(df_return[df_return['Gene'] == i]['Z-Score'].values.tolist()[0]) if i in df_return['Gene'].values.tolist() else 0 for i in clust2]
            n_colors_dict = dict(zip(clust2, n_colors))
        elif color_by == 'degree_centrality':
            n_colors_dict = nx.degree_centrality(subgraph)
        elif color_by == 'closeness_centrality':
            n_colors_dict = nx.closeness_centrality(subgraph)
        elif color_by == 'betweenness_centrality':
            n_colors_dict = nx.betweenness_centrality(subgraph)
        
        if color_by == 'ceres_score' or color_by == 'z_score':
            
            test_c = np.linspace(min(n_colors)-0.1, max(n_colors)+0.1, cres)
            cmaph = plt.cm.get_cmap(cmap, cres)
            cmaph = [(255*cmaph(i)[0], 255*cmaph(i)[1], 255*cmaph(i)[2])  for i in range(cres)]

            def color_this(fl, breakpoints=test_c, cat=cmaph):
                return cat[bisect(breakpoints, fl)]
                    
            nodes = [
                {
                    'data': {'id': label, 'label': label, 'size': sizes_dict[label], 'color': color_this(n_colors_dict[label])},
                    'classes': n_class_dict[label],
                }
                for label in clust2
            ]
            
        elif color_by == 'degree_centrality' or color_by == 'betweenness_centrality' or color_by == 'closeness_centrality':
            
            low = min([v for k, v in n_colors_dict.items()])
            high = max([v for k, v in n_colors_dict.items()])
            percent_10 = high - low / 10
            test_c = np.linspace(low-percent_10, high+percent_10, cres)
            cmaph = plt.cm.get_cmap(cmap, cres)
            cmaph = [(255*cmaph(i)[0], 255*cmaph(i)[1], 255*cmaph(i)[2])  for i in range(cres)]

            def color_this(fl, breakpoints=test_c, cat=cmaph):
                return cat[bisect(breakpoints, fl)]
                    
            nodes = [
                {
                    'data': {'id': label, 'label': label, 'size': sizes_dict[label], 'color': color_this(n_colors_dict[label])},
                    'classes': n_class_dict[label],
                }
                for label in clust2
            ]
            
        elif color_by == 'community':
            
            nodes = [
                {
                    'data': {'id': label, 'label': label, 'size': sizes_dict[label], 'color': n_colors_dict[label]},
                    'classes': n_class_dict[label],
                }
                for label in clust2
            ]

        edges = [
            {'data': {'source': source, 'target': target},
             'classes': e_class_dict[(source, target)]
            }
            for source, target in G.edges(clust)
        ]

        elements = nodes + edges
        
        if size_by == 'ceres_score':
            group_select_dict = {
                            'selector': 'node',
                            'style': {
                                'content': 'data(label)',
                                'width': 'mapData(size,  -3, 1, 100, 20)',
                                'height': 'mapData(size,  -3, 1, 100, 20)',
                                'background-color': 'data(color)',
                            }
                        }
                
        elif size_by == 'z_score':
            group_select_dict = {
                            'selector': 'node',
                            'style': {
                                'content': 'data(label)',
                                'width': 'mapData(size,  0, 2, 20, 100)',
                                'height': 'mapData(size,  0, 2, 20, 100)',
                                'background-color': 'data(color)',
                            }
                        }
                
        elif size_by == 'degree_centrality' or size_by == 'betweenness_centrality' or size_by == 'closeness_centrality':
            low_s = min([v for k, v in sizes_dict.items()])
            high_s = max([v for k, v in sizes_dict.items()])
            map_data = 'mapData(size, ' + str(low_s) + ', ' + str(high_s) + ', 20, 100)'
            group_select_dict = {
                            'selector': 'node',
                            'style': {
                                'content': 'data(label)',
                                'width': map_data,
                                'height': map_data,
                                'background-color': 'data(color)',
                            }
                        }

        network =  [html.Br(),
                    dbc.Row([
                        dbc.Col(width=2),
                        dbc.Col(dbc.Label("Network Topology of Selected Communities", style={'text-align': 'center', 'font-size': '150%', "font-weight": "bold", 'font-family': 'Arial'}), width=9),
                        dbc.Col(width=1),
                    ]),
                    html.Br(),
                    cyto.Cytoscape(
                    id='cytoscape-layout-2',
                    elements=elements,
                    style={'width': '100%', 'height': '500px'},
                    layout={
                        'name': layout,
                        'animate': True
                    },
                    stylesheet=[
                        # Group selectors
                        group_select_dict,

                        # Class selectors
                        #{
                        #    'selector': '.red',
                        #    'style': {
                        #        'background-color': 'magenta',
                        #        'line-color': 'red'
                        #    }
                        #},
                        #{
                        #    'selector': '.blue',
                        #    'style': {
                        #        'background-color': 'lightblue',
                        #        'line-color': 'grey'
                        #    }
                        #}
                    ]
                )]

        df_return['Gene'] = df_return['Gene'].apply(lambda x: create_gene_link(x, 'uniprot', 'gene'))

        return network, dbc.Table.from_dataframe(df_return, striped=True, bordered=True, hover=True, size='sm')

@app.callback(Output('cluster_graph_1', 'figure'),
              Output('cluster_feats_1', 'figure'),
              Output('select_cluster_2', 'options'),
              Input('select_community_2', 'value'),
              Input('select_cluster_layout', 'value'),
              Input('select_node_colormap_2', 'value'))
def create_cluster_graph(landmark, layout, cmap):
    if landmark is None or layout is None or cmap is None:
        return dash.no_update
    elif landmark is not None and layout is not None:

        # for s3:
        #phate_op_key = "cluster_models/" + 'louvain_' + str(landmark) + "_phate_op.pkl"
        #phate_op = pickle.load(return_s3_object(s3, s3_bucket_name, phate_op_key))
        
        #coords_key = "cluster_models/" + 'louvain_' + str(landmark) + "_coords.npy"
        #coords = np.load(return_s3_object(s3, s3_bucket_name, coords_key))

        #clusters = phate.cluster.kmeans(phate_op, k=7)
        
        #cluster_cell_lines_key = "cluster_models/" + 'louvain_' + str(landmark) + '_cell_lines.pkl'
        #cluster_cell_lines = pickle.load(return_s3_object(s3, s3_bucket_name, cluster_cell_lines_key))

        # for local:
        phate_op_key = "./ceres-infer/cluster_models/" + 'louvain_' + str(landmark) + "_phate_op.pkl"
        with open(phate_op_key, 'rb') as f:
            phate_op = pickle.load(f)
        
        coords_key = "./ceres-infer/cluster_models/" + 'louvain_' + str(landmark) + "_coords.npy"
        coords = np.load(coords_key)

        clusters = phate.cluster.kmeans(phate_op, k=7)
        
        cluster_cell_lines_key = "./ceres-infer/cluster_models/" + 'louvain_' + str(landmark) + '_cell_lines.pkl'
        with open(cluster_cell_lines_key, 'rb') as f:
            cluster_cell_lines = pickle.load(f)

        Y_phate_tsne = coords[int(layout[0])]

        fig = go.Figure()

        unique = set(list(clusters))
        cmaph = plt.cm.get_cmap(cmap, len(unique))
        color_bys = unique
        
        cmaph = [(cmaph(i)[0], cmaph(i)[1], cmaph(i)[2], cmaph(i)[3])  for i in range(len(unique))]
        cmaph = [mpl.colors.rgb2hex(i, keep_alpha=False) for i in cmaph]

        color_by_dict = dict(zip(color_bys, cmaph))

        for clust in unique:
            fig.add_trace(go.Scatter(x=[row[0] for row in Y_phate_tsne[clusters==clust]], 
                                    y=[row[1] for row in Y_phate_tsne[clusters==clust]],
                                    mode='markers',
                                    marker=dict(
                                        color=color_by_dict[clust],
                                        size=10,
                                        line=dict(
                                        color='DarkSlateGray',
                                        width=1
                                        )),
                                    text=[cluster_cell_lines[i] for i in range(len(cluster_cell_lines)-1) if [clusters==clust][0][i]==True],
                                    hoverinfo='text',
                                    name = "cluster " + str(clust)))
        
        fig.update_yaxes(zeroline=True, 
                         zerolinecolor='rgb(211,211,211,0.8)', 
                         zerolinewidth=0.7, 
                         showgrid=False, 
                         title="<b>" + layout[1:] + "2",
                         tickfont_family="Arial Black", 
                         tickfont=dict(size=18),
                         showticklabels=False)
        
        fig.update_xaxes(zeroline=True, 
                         zerolinecolor='rgb(211,211,211,0.8)', 
                         zerolinewidth=0.7, 
                         showgrid=False, 
                         title="<b>" + layout[1:] + "1", 
                         tickfont_family="Arial Black", 
                         tickfont=dict(size=18),
                         showticklabels=True)
        
        fig.update_layout(
        plot_bgcolor='white',
        title_text = '<b>PHATE Clusters for Louvain Community: ' + str(landmark) + '</b> <br>' + layout[1:] + ' Projection',
        title_x = 0.5,
        font_family = 'Arial',
        font_color = 'black',
        font_size=16,
        legend = dict(font=dict(size=16, family="Arial", color="black"), orientation='h', yanchor='top', y=-0.1),
        margin=dict(t=80, b=80),
        height=700)

        # for s3:
        #feature_key = "cluster_models/" + 'louvain_' + str(landmark) + "_feat_importance_rf.pkl"
        #feature_dict = pickle.load(return_s3_object(s3, s3_bucket_name, feature_key))
        
        # for local:
        feature_key = "./ceres-infer/cluster_models/" + 'louvain_' + str(landmark) + "_feat_importance_rf.pkl"
        with open(feature_key, 'rb') as f:
            feature_dict = pickle.load(f)

        fig2 = make_subplots(rows=2, cols=4, x_title='<b>Feature Score')

        unique = ['0', '1', '2', '3', '4', '5', '6']
        cmaph = plt.cm.get_cmap(cmap, len(unique))
        
        cmaph = [(cmaph(i)[0], cmaph(i)[1], cmaph(i)[2], cmaph(i)[3])  for i in range(len(unique))]
        cmaph = [mpl.colors.rgb2hex(i, keep_alpha=False) for i in cmaph]

        color_by_dict = dict(zip(unique, cmaph))
        
        clust_coord_dict = {'0': (1, 1), '1': (1, 2), '2': (1, 3), '3': (1, 4), '4': (2, 1), '5': (2, 2), '6': (2, 3)}

        for clust in unique:
            importance_dict = feature_dict[int(clust)]
            feat_labels = list(importance_dict.keys())
            importances = list(importance_dict.values())
            importances, feat_labels = zip(*sorted(zip(importances, feat_labels), reverse=True))
            fig2.add_trace(go.Bar(
                y=feat_labels,
                x=importances,
                name='cluster ' + clust,
                orientation='h',
                marker=dict(
                    color=color_by_dict[clust],
                    line=dict(
                        color='black',
                        width=0.2
                        ),)
                ), row = clust_coord_dict[clust][0], col = clust_coord_dict[clust][1])
            fig2.update_yaxes(zeroline=True, zerolinecolor='black', zerolinewidth=2, 
                        tickfont_family="Arial Black", tickfont=dict(size=12, color='black'),showgrid=False, showticklabels=True)
            fig2.update_xaxes(zeroline=True, zerolinecolor='black', zerolinewidth=2, 
                        tickfont_family="Arial Black", tickfont=dict(size=12, color='black'),showgrid=False, showticklabels=True)

        fig2.update_layout(plot_bgcolor='white',
                            height=700,
                            title_text = '<b>Feature Importance Scores for Clusters of Louvain Community: ' + str(landmark),
                            title_x = 0.5,
                            font_family = 'Arial',
                            font_color = 'black',
                            font_size=16,
                            legend = dict(font=dict(size=16, family="Arial", color="black")))
        fig2.layout.annotations[0]["font"] = {'size': 20, 'family': 'Arial', 'color':'black'}
        
        return fig, fig2, [{'label': "cluster " + str(value), 'value': value} for value in unique]

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
        # for s3:
        #feature_key = "cluster_models/" + 'louvain_' + str(landmark) + "_feat_importance_rf.pkl"
        #feature_dict = pickle.load(return_s3_object(s3, s3_bucket_name, feature_key))
        
        # for local:
        feature_key = "./ceres-infer/cluster_models/" + 'louvain_' + str(landmark) + "_feat_importance_rf.pkl"
        with open(feature_key, 'rb') as f:
            feature_dict = pickle.load(f)
        
        importance_dict = feature_dict[int(clust)]
        feat_labels = list(importance_dict.keys())
        importances = list(importance_dict.values())
        importances, feat_labels = zip(*sorted(zip(importances, feat_labels), reverse=True))

        cluster = modularity[landmark]
        
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
        
        return dbc.Table.from_dataframe(df_return, striped=True, bordered=True, hover=True, size='sm')

app.register_celery_tasks()

if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8080, debug=True)
