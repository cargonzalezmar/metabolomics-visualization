#!/usr/bin/env python

import plotly.graph_objects as go
import plotly.express as px


template = "simple_white"


##############################################################################
#############                 Consensus Graphs           #####################
##############################################################################


def init_consensus_graph(df):   
    fig = px.scatter(
                    df, 
                    x="RT", 
                    y="mz",
                    marginal_x="violin", 
                    marginal_y="violin",
                    # height=600, width=800,
                    template=template
                )

    fig.update_layout(
                    showlegend=False,
                    #title_text=f"RT vs MZ", 
                    title_x=0.5
                )
    return fig

def create_consensus_graph(df1, df2, sample_name):

    fig = px.scatter(
                df1, 
                x="RT", 
                y="mz",
                marginal_x="violin", 
                marginal_y="violin", 
                color=sample_name,
                # height=600, width=800,
                template=template
            )

    fig.add_trace(
                go.Scatter(
                x=df2["RT"],
                y=df2["mz"],
                mode='markers',
                marker=dict(color="rgb(255,255,255)", line=dict(width=1, color='Grey'))
                
            ))

    fig.update_layout(
                    showlegend=False,
                    # title_text=f"RT vs MZ by intensity on {sample_name}", 
                    title_x=0.5
                )

    fig.layout.coloraxis.colorbar.title = 'Intensity'

    return fig


##############################################################################
#############                   MS1 Graphs               #####################
##############################################################################


def create_feature_graph(feat_maps_df):

    fig = px.line(
            feat_maps_df, 
            x="mz", 
            y="intensity",
            height=500, width=900,
            template = template,
            hover_name = feat_maps_df.index
        )

    fig.update_layout(
            showlegend=False,
            title_text=f"mz vs intensity", 
            title_x=0.5
        )

    fig.update_traces(customdata=feat_maps_df.index)

    return fig 


##############################################################################
#############                   MS2 Graphs               #####################
##############################################################################


def create_ms2_graph(df, index_list, sample_name):

    data_fig = []
    for index in index_list:
        rt = df[sample_name][index]["RT"]
        mz_int_data = df[sample_name][index]["mz_int_data"]

        trace = go.Line(
            x = mz_int_data["mz"],
            y = mz_int_data["intensity"],
            name = rt
        )
        data_fig.append(trace)

    fig = go.Figure(data_fig)
    
    fig.update_layout(
        template=template, 
        title="MS2",
        legend_title="Retention Time",
        title_text= f"MZ vs intensity {sample_name}", 
        title_x=0.5,
        height=500, width=900,
        )

    return 

