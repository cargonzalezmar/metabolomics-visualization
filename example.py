import os
import pandas as pd
import numpy as np

import plotly.graph_objects as go
import plotly.express as px

from modules.graphs import *

from dash import Dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import dash_bio as dashbio

results_directory = "subset_results"

consensus = pd.read_pickle(os.path.join(results_directory, "FeatureMatrix_formulas_structures.pkl"))
samples = [col[:-5] for col in consensus.columns if col.endswith("mzML")]
consensus.columns = [col.split(".mzML")[0] if ".mzML" in col else col for col in consensus.columns]
spectra = {sample: pd.read_pickle(os.path.join(results_directory, "ms_df", sample+".pkl")) for sample in samples}
features = {sample: pd.read_pickle(os.path.join(results_directory, "feature_map_df", sample+"fm.pkl")) for sample in samples}

sample_dict = [dict(label=x, value=x) for x in samples]
sample_dict.append({'label':'All Samples', 'value':'All Samples'})

app = Dash(__name__)

app.layout = html.Div([html.Div([html.Div([html.H1(children='Metabolomic Visualization', id='title'),
                                           html.P(children="Holiiiii Axeeeeel!!!!")],
                                           id='left-header'), 
                                 html.A([html.Img(src='./assets/git.png', id='logo'),
                                         html.P(children="GitHub")],
                                        href='https://github.com/cargonzalezmar/metabolomics-visualization',
                                        id='right-header')],
                                        id='header'),
                      html.Div([dcc.Dropdown(options=sample_dict, id='samples_dropdown', value='All Samples'),
                                dcc.Graph(id='consensus_graph'),
                                dcc.Graph(id='feature_spectra'),
                                dcc.Graph(figure=dashbio.Clustergram(
                                    data=consensus[samples].loc[list(consensus[samples].index)].values,
                                    column_labels=samples,
                                    row_labels=list(consensus[samples].index),
                                    color_threshold={
                                    'row': 250,
                                    'col': 700
                                    },
                                    hidden_labels='row',
                                    height=800,
                                    width=700
                                    ))
                                ], 
                                id='body'),
                      html.Div(children='Axel Walter & Carolina González',id='footer')  
                                ]
                     )
                     
@app.callback(
    Output("consensus_graph", "figure"), 
    [Input("samples_dropdown", "value")])
def update_consensus_graph(sample_name):
    if sample_name=='All Samples':
        fig = init_consensus_graph(consensus)
    else:
        df1 = consensus[consensus[sample_name]>0]
        df2 = consensus[consensus[sample_name]==0]
        fig = create_consensus_graph(df1, df2, sample_name)
    return fig

@app.callback(
    Output("feature_spectra", "figure"), 
    [Input("samples_dropdown", "value")])
def update_ms1_graph(sample_name): 
    feat_maps_df = features[sample_name].sort_values(by='mz')
    fig = create_feature_graph(feat_maps_df)
    return fig

if __name__=='__main__':
    app.run_server(debug=True)
