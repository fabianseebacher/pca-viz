import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import numpy as np
import pandas as pd
import yaml

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

df = pd.read_csv("N:\\IDO_Proteomics_CellBiol\\Fabian\\PCA_Classification_f6.txt", sep="\t", skiprows=[1,2,3])
# df["Winner"] = df["Winner"].replace(np.nan, "None")
# df['Organell_f6'] = df['Organell_f6'].replace(np.nan, "None")

### LAYOUT ###
    
app.layout = html.Div([

    # html.H1('PCA Vizualizer'),

    html.Div([
        
        html.Div([
            html.P("Enter path"),
            dcc.Input(
                id="pathtext",
                type="text",
                value=r"N:\\IDO_Proteomics_CellBiol\\Fabian\\PCA_Classification_f6.txt"
            )

        ])    
    ]),

    html.Details([
            html.Summary('Options'),
            html.Div([
                html.Div([
                    html.P("Select first component:"),
                    dcc.Dropdown(
                        id='pca-1',
                        options=[{'label': i, 'value': i} for i in df.columns],
                        value='Component 1'
                    )
                ], style={'display': 'inline-block', 'width': "33%"}),
                html.Div([
                    html.P("Select second component:"),
                    dcc.Dropdown(
                        id='pca-2',
                        options=[{'label': i, 'value': i} for i in df.columns],
                        value='Component 2'
                    )
                ], style={'display': 'inline-block', 'width': "33%"}),
                html.Div([
                    html.P("Select third component:"),
                    dcc.Dropdown(
                        id='pca-3',
                        options=[{'label': i, 'value': i} for i in df.columns],
                        value='Component 3'
                    )
                ], style={'display': 'inline-block', 'width': "33%"}),
            ]),
            html.Div([
                html.P("Select column with marker organelle annotation:"),
                dcc.Dropdown(
                    id='marker-col',
                    options=[{'label': i, 'value': i} for i in df.columns],
                    value='Organell_f6'
                )
            ], style={'display': 'inline-block', 'width': "48%"}),

            html.Div([
                html.P("Select column with predicted organelle :"),
                dcc.Dropdown(
                    id='predict-col',
                    options=[{'label': i, 'value': i} for i in df.columns],
                    value='Winner'
                )
            ], style={'display': 'inline-block', 'width': "48%"})
    ]),
        
    html.Div([
        html.P("Select protein:"),
        dcc.Dropdown(
            id='proteinselect',
            options=[{'label': i, 'value': i} for i in df["Genes"].sort_values()],
            value=None,
            multi=True
        )
    ], style={'display': 'inline-block', 'width': "48%"}),

    html.Div([
        html.P("Color mode"),
        dcc.Dropdown(
            id='colormode',
            options=[{'label': i, 'value': i} for i in ['None', 'markers only', 'all predicted']],
            value='all predicted'
        )
    ], style={'display': 'inline-block', 'width': "48%"}),

    html.Br(),


    html.Div([
        dcc.Graph(id='3d_pca'),
    ]),

    html.Label(id='numbers', children='Display numbers here')
])

### CALLBACKS ###

@app.callback(
    Output('3d_pca', 'figure'), Output('numbers', 'children'),
    Input('pathtext', 'value'), Input('colormode', 'value'), Input('proteinselect', 'value'),
    Input('marker-col', 'value'), Input('predict-col', 'value'), Input('pca-1', 'value'), Input('pca-2', 'value'), Input('pca-3', 'value'))
def update_graph(pathtext, colormode, proteins, mcol, pcol, pc1, pc2, pc3):
    print()
    df = pd.read_csv(pathtext, sep="\t", skiprows=[1,2,3])
    
    df[pcol] = df[pcol].replace(np.nan, "None")
    df[mcol] = df[mcol].replace(np.nan, "None")
    
    if colormode == 'all predicted':
        fig = px.scatter_3d(df, x=pc1, y=pc2, z=pc3, template="plotly_white", color=pcol, opacity=0.5, hover_name="Genes", color_discrete_sequence=["cyan", "lightgrey", "orange", "purple", "yellow", "red", "green", "blue", "lightgreen", 'magenta', 'black'])
    elif colormode == 'markers only':
        fig = px.scatter_3d(df, x=pc1, y=pc2, z=pc3, template="plotly_white", color=mcol, opacity=0.5, hover_name="Genes", color_discrete_sequence=["lightgrey", "orange", "red", "green", "cyan", "purple", "blue", "lightgreen", "yellow", 'magenta', 'black' ]) 
    else:
        fig = px.scatter_3d(df, x=pc1, y=pc2, z=pc3, template="plotly_white", opacity=0.5, hover_name="Genes")
        fig.update_traces(marker=dict(color="lightgrey"),
                  selector=dict(mode='markers'))
    
    
    
    fig.update_traces(marker=dict(size=3),
              selector=dict(mode='markers'))

    if proteins:
        df_subset = df[df['Genes'].isin(proteins)]
        fig2 = px.scatter_3d(df_subset, x=pc1, y=pc2, z=pc3, template="plotly_white", hover_name="Genes")
        fig2.update_traces(marker=dict(size=6, color="red"),
                  selector=dict(mode='markers'))
        fig.add_trace(fig2.data[0])
    
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0), legend=dict(
        itemsizing = 'constant',
        orientation = 'h',
        yanchor = 'bottom',
        y = -0.02,
        xanchor = 'right',
        x = 1 
     ) ,uirevision=True)
    
    stats = []
    #stats.append("Markers:")
    for i in df[mcol].unique():
        stats.append(str(i) + ": " + str( len(df[df[mcol] == i]) ))
    
    stats2 = []
    #stats2.append("Predicted:")
    for i in df[pcol].unique():
        stats2.append(str(i) + ": " + str( len(df[df[pcol] == i]) ))
    return [fig, [html.B("Markers: "), ' '.join(map(str, stats)), html.Br(), html.B("Predicted: "),' '.join(map(str, stats2))]]

# Update dropdowns for available columns based on path provided
@app.callback(
    Output('marker-col', 'options'), Output('predict-col', 'options'), Output('proteinselect', 'options'), Output('pca-1', 'options'), Output('pca-2', 'options'), Output('pca-3', 'options'),
    Input('pathtext', 'value'))
def set_dir_options(selected_path):
    df = pd.read_csv(selected_path, sep="\t", skiprows=[1,2,3])
    opt = [{'label': i, 'value': i} for i in df.columns]
    prot_opt = [{'label': i, 'value': i} for i in df["Genes"].sort_values()]
    return [ opt, opt, prot_opt, opt, opt, opt ]



if __name__ == '__main__':
    with open("settings.yml", 'r') as stream:
        settings = yaml.safe_load(stream)
    print(settings)
    app.run_server(debug=False, host=settings['ip'], port=int(settings['port']) )