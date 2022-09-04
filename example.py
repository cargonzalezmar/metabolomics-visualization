import os
import pandas as pd
import numpy as np

import plotly.graph_objects as go
import plotly.express as px

from dash import Dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output

results_directory = "subset_results"

consensus = pd.read_pickle(os.path.join(results_directory, "FeatureMatrix_formulas_structures.pkl"))
samples = [col[:-5] for col in consensus.columns if col.endswith("mzML")]
spectra = {sample: pd.read_pickle(os.path.join(results_directory, "ms_df", sample+".pkl")) for sample in samples}
features = {sample: pd.read_pickle(os.path.join(results_directory, "feature_map_df", sample+"fm.pkl")) for sample in samples}

sample_dict = [dict(label=x, value=x) for x in samples]

app = Dash(__name__)

app.layout = html.Div([html.Div([html.Div([html.H1(children='Metabolomic Visualization', id='title'),
                                           html.P(children="Holiiiii Axeeeeel!!!!")],
                                           id='left-header'), 
                                 html.A([html.Img(src='./assets/git.png', id='logo'),
                                         html.P(children="GitHub")],
                                        href='https://github.com/cargonzalezmar/metabolomics-visualization',
                                        id='right-header')],
                                 id='header'),
                      html.Div([dcc.Dropdown(options=sample_dict, id='samples_dropdown', value=samples[0]),
                                dcc.Graph(id='ms1_spectra'),
                                ], 
                                id='body'),
                      html.Div(children='Axel Walter & Carolina González',id='footer')  
                                ]
                     )

@app.callback(
    Output("ms1_spectra", "figure"), 
    [Input("samples_dropdown", "value")])
def update_ms1_graph(samples):
    features_df = features[samples].sort_values(by='mz')
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=features_df['mz'],
                     y=features_df['intensity'],
                     line={'color':'black', 'width':0.5}
                    ))
    fig.update_layout(
            showlegend=False,
            title_text=f"mz vs intensity", 
            title_x=0.5
        )
    fig.update_traces(customdata=features_df.index)
    return fig 
    # fig = px.line(
    #         features[samples].sort_values(by='mz'), 
    #         x="mz", 
    #         y="intensity",
    #         height=500, width=900,
    #         template = 'simple_white'
    #     )

    # fig.update_layout(
    #         showlegend=False,
    #         title_text=f"mz vs intensity", 
    #         title_x=0.5
    #     )

    # return fig 


if __name__=='__main__':
    app.run_server(debug=True)
features_df