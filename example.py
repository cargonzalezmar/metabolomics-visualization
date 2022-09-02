from dash import Dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output

app = Dash(__name__)

app.layout = html.Div(
                [html.Div([html.H1(children='Metabolomic Visualization', id='title')],
                id='header')]
            )

if __name__=='__main__':
    app.run_server(debug=True)
