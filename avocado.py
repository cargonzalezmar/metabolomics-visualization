import os

from dash import Dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output

import pandas as pd

results_directory = "results"

# load all the data into consensus (dataframe), spectra and features (dict sample names as keys, dataframes as values)
consensus = pd.read_feather(os.path.join(results_directory, "FeatureMatrix.ftr"))
samples = [col[:-5] for col in consensus.columns if col.endswith("mzML")]
spectra = {sample: pd.read_feather(os.path.join(results_directory, "ms_df", sample+".ftr")) for sample in samples}
features = {sample: pd.read_feather(os.path.join(results_directory, "feature_map_df", sample+".ftr")) for sample in samples}

# initialize the app
app = Dash(__name__)

# define app layout
app.layout = html.Div(
    children=[
        dcc.Dropdown(
            options=samples+["all"],
            id="consensus-sample",
        ),
        dcc.Graph("consensus-plot")   
    ]
)

# draw consensus map graph
@app.callback(Output("consensus-plot", "figure"), [Input("consensus-sample", "value")])
def update(sample):
    figure={
        "data": [
            {
                "x": [spectra[sample][i]["RT"] for i in data.dct_ms_df[sample].keys()],
                "y": [max(spectra[sample][i]['mz_int_data']["intensity"]) for i in spectra[sample].keys()],
                "type": "line",
                "mode": "lines",
                "marker": {"color": "rgba(0, 0, 0, 0.8)"}
            }
        ],
        "layout": {
            "title": sample,
            "xaxis": {"title": "mz"},
            "yaxis": {"title": "intensity"},
            "showlegend": False
        }
    }
    return figure

# draw consensus hovered MS2 graph

if __name__ == "__main__":
    app.run_server(debug=True)