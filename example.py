from crypt import METHOD_SHA512
import os
from turtle import title
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

def plot_ms_spectrum(df_spectrum, title):
    """
    Takes a pandas Dataframe with one row (spectrum) and generates a needle plot with m/z and intensity dimension.
    """
    def create_spectra(x,y, zero=0):
        x=np.repeat(x,3)
        y=np.repeat(y,3)
        y[::3]=y[2::3]=zero
        return x,y

    x,y = create_spectra(df_spectrum['mzarray'].tolist()[0], df_spectrum['intarray'].tolist()[0])
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x,y=y))
    fig.update_traces(showlegend=True)
    fig.update_layout(
        showlegend=False,
        title_text=title,
        xaxis_title="m/z",
        yaxis_title="intensity"
    )
    return fig

app = Dash(__name__)

app.layout = html.Div([html.Div([html.Div([html.H1(children='Metabolomic Visualization', id='title'),
                                           html.P(children=" Holiiiii Axeeeeel!!!!")],
                                           id='left-header'), 
                                 html.A([html.Img(src='./assets/git.png', id='logo'),
                                         html.P(children="GitHub")],
                                        href='https://github.com/cargonzalezmar/metabolomics-visualization',
                                        id='right-header')], id='header'),
                      html.Div([html.Br(),
                                dcc.Dropdown(options=sample_dict, id='samples_dropdown', value='All Samples'),
                                html.Br(),
                                html.H2(children="Consensus features", className="subtitle"),
                                html.Div([
                                    dcc.Graph(id='consensus_graph', className="panel-item-left-60"),
                                    dcc.Graph(figure=dashbio.Clustergram(
                                        data=consensus[samples].loc[list(consensus[samples].index)].values,
                                        column_labels=samples,
                                        row_labels=list(consensus[samples].index),
                                        color_map="viridis",
                                        color_threshold={
                                        'row': 250,
                                        'col': 700
                                        },
                                        hidden_labels='row',
                                        ), className="panel-item-right-40")
                                ], className="dual-panel-60-40"),
                                html.H2(children="Features", className="subtitle"),
                                dcc.Graph(id="feature-scatter"),
                                html.Div([
                                    dcc.Graph(id="feature-chromatogram", className="panel-item-left"),
                                    dcc.Graph(id="feature-ms2", className="panel-item-right")
                                ], className="dual-panel"),
                                html.H2(children="MS1 data", className="subtitle"),
                                html.Div([dcc.Graph(id="BPC", className="panel-item-left"),
                                        dcc.Graph(id="ms1-spectra", className="panel-item-right")], id='ms1-panel', className="dual-panel"),
                                html.H2(children="MS2 data", className="subtitle"),
                                html.Div([dcc.Graph(id="ms2-scatter-plot", className="panel-item-left"), 
                                        dcc.Graph(id="ms2-spectrum", className="panel-item-right")], id='ms2-panel', className="dual-panel")
                                ], 
                                id='body'),
                      html.Div(children='Axel Walter & Carolina González',id='footer')  
                                ]
                     )
                     
@app.callback(
    Output("consensus_graph", "figure"), 
    [Input("samples_dropdown", "value")])
def update(sample_name):
    if sample_name=='All Samples':
        fig = init_consensus_graph(consensus)
    else:
        df1 = consensus[consensus[sample_name]>0]
        df2 = consensus[consensus[sample_name]==0]
        fig = create_consensus_graph(df1, df2, sample_name)
    return fig

@app.callback(
    Output("feature-scatter", "figure"),
    [Input("samples_dropdown", "value")]
)
def update(sample):
    if sample=="All Samples":
        return go.Figure()
    fig = px.scatter(features[sample], x="RT", y="mz", title=f"{sample} features", color="intensity", color_continuous_scale="Hot_r")
    # fig.update_traces(marker=dict(color="purple"))
    fig.update_layout(yaxis={"title": "m/z"})
    return fig

@app.callback(
    Output("feature-chromatogram", "figure"),
    [Input("samples_dropdown", "value"), Input("feature-scatter", "clickData")])
def update(sample, clickData):
    if sample=="All Samples":
        return go.Figure()
    if clickData:
        df = features[sample].loc[(features[sample]["RT"] == clickData["points"][0]["x"]) & 
                                (features[sample]["mz"] == clickData["points"][0]["y"])]
        title = "feature chromatogram @"+ str(round(clickData["points"][0]["x"])) + " RT @" + str(round(clickData["points"][0]["y"], 5)) +" m/z"
        fig = px.line(x=df["chroms"][0][0]["rt"], y=df["chroms"][0][0]["intensity"], title=title)
        # add second isotope trace chromatogram if exists
        if len(df["chroms"][0]) > 1:
            fig2 = px.line(x=df["chroms"][0][1]["rt"], y=df["chroms"][0][1]["intensity"])
            fig.add_trace(fig2.data[0])
        fig.update_layout(xaxis={"title": "RT"}, yaxis={"title": "intensity"})
        return fig
    else:
        return go.Figure()

@app.callback(
    Output("feature-ms2", "figure"),
    [Input("samples_dropdown", "value"), Input("feature-scatter", "clickData")])
def update(sample, clickData):
    if sample=="All Samples":
        return go.Figure()
    if clickData:
        ms2_spectra = features[sample].loc[(features[sample]["RT"] == clickData["points"][0]["x"]) & 
                                        (features[sample]["mz"] == clickData["points"][0]["y"])]["MS2_spectra"][0]
        fig = go.Figure()
        while len(ms2_spectra) > 0:
            fig.add_trace(plot_ms_spectrum(spectra[sample].iloc[ms2_spectra.pop()].to_frame().transpose(), title="MS2 spectrum").data[0])
        title = "annotated MS2 spectra @"+ str(round(clickData["points"][0]["x"])) + " RT @" + str(round(clickData["points"][0]["y"], 5)) +" m/z"
        fig.update_layout(title=title)
        return fig
    else:
        return go.Figure()

# draw BPC of selected sample
@app.callback(Output("BPC", "figure"), [Input("samples_dropdown", "value")])
def update(sample):
    if sample=="All Samples":
        return go.Figure()
    if sample == "all":
        return {}
    df_ms1 = spectra[sample].loc[spectra[sample]["mslevel"] == 1]
    fig = px.line(df_ms1, x="RT", y=[max(intensity_array) for intensity_array in df_ms1["intarray"]], title=f"{sample} BPC", labels={"y": "intensity (cps)"})
    fig.data[0].update(mode='markers+lines', marker_symbol="x", marker=dict(color="orange"))
    return fig

@app.callback(Output("ms1-spectra", "figure"), [Input("samples_dropdown", "value"), Input("BPC", "clickData")])
def update(sample, clickData):
    if sample=="All Samples":
        return go.Figure()
    if clickData:
        df_ms1 = spectra[sample].loc[spectra[sample]["mslevel"] == 1]
        df_spectrum = df_ms1.loc[df_ms1["RT"] == clickData["points"][0]["x"]]
        return plot_ms_spectrum(df_spectrum, sample + " MS1 @" + str(round(clickData["points"][0]["x"])) + " s")
    else:
        return go.Figure()

# draw scatter plot of MS2 precursors
@app.callback(Output("ms2-scatter-plot", "figure"), [Input("samples_dropdown", "value")])
def update(sample):
    if sample=="All Samples":
        return go.Figure()
    df_ms2 = spectra[sample].loc[spectra[sample]["mslevel"] == 2]
    fig = px.scatter(df_ms2, x="RT", y="precursors", title=f"{sample} MS2 precursors", labels={"y": "precursor m/z"})
    fig.update_traces(marker=dict(color="orange"))
    fig.update_layout(yaxis={"title": "precursor m/z"})
    return fig

# draw MS2 spectrum on hover of ms2-scatter-plot
@app.callback(Output("ms2-spectrum", "figure"), [Input("samples_dropdown", "value"), Input("ms2-scatter-plot", "clickData")])
def update(sample, clickData):
    if sample=="All Samples":
        return go.Figure()
    if clickData:
        df_ms2 = spectra[sample].loc[spectra[sample]["mslevel"] == 2]
        df_spectrum = df_ms2.loc[df_ms2["RT"] == clickData["points"][0]["x"]]
        return plot_ms_spectrum(df_spectrum, sample + " MS2 @" + str(round(clickData["points"][0]["x"])) + " s @" + str(round(clickData["points"][0]["y"], 5)) + " precursor m/z")
    else:
        return go.Figure()

if __name__=='__main__':
    app.run_server(debug=True)
