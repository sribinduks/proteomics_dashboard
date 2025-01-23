import pandas as pd
import plotly.express as px
import numpy as np
import dash
from dash import dcc
from dash import html

# Create dummy data
np.random.seed(42)
df = pd.DataFrame({
    'uniprot': ['P00001', 'P00002', 'P00003', 'P00004', 'P00005'],
    'logfc': np.random.normal(0, 1, 5),
    '-logp': np.abs(np.random.normal(0, 1, 5))
})

# Threshold values
threshold_logfc = 1.5
threshold_p = 1

# Dictionary relating some uniprots  
relations = {'P00001': ['P00002'],
             'P00004': ['P00005']}

app = dash.Dash(__name__)

fig = px.scatter(df, x='logfc', y='-logp', color_discrete_sequence=['grey']*len(df))

app.layout = html.Div([
    dcc.Graph(id='graph', figure=fig),
    html.Div(id='click-data') 
])

@app.callback(
    dash.dependencies.Output('graph', 'figure'),
    [dash.dependencies.Input('graph', 'clickData'),
     dash.dependencies.State('graph', 'figure')])
def update_data(clickData, figure):
    changed_points = []
    if clickData:
        sel_uniprot = df.loc[clickData['points'][0]['pointIndex'], 'uniprot']
        changed_points.append(clickData['points'][0]['pointIndex'])
        for i, row in df.iterrows():
            if row['uniprot'] == sel_uniprot or row['uniprot'] in relations.get(sel_uniprot, []):
                changed_points.append(i)

    colors = ['grey']*len(df)
    for i in changed_points:
        colors[i] = 'green'
        
    opacity = [0.25]*len(df)
    for i in changed_points:
        opacity[i] = 1.0
        
    df['colors'] = colors
    df['opacity'] = opacity
    
    figure['data'][0].data= df
    figure['data'][0].marker.color = df['colors']
    figure['data'][0].marker.opacity = df['opacity']

    return figure

if __name__ == '__main__':
    app.run_server(debug=True)