# File: app_sri.py (newest)

import dash
import math
from dash import dcc
from dash import html
from dash import State
from dash.dependencies import Input, Output
import pandas as pd
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import plotly.express as px
import json
import numpy as np
import plotly.graph_objects as go
from data_loader import load_data, load_txt
from get_molecule_smile import get_smile_url

# external stylesheet
external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

# Load the data
folder = 'data'
data = load_data(folder)
txt_data = load_txt(folder)
hits = data["pde3a_inh_ip_hits_p0001.csv"]

network_df = pd.read_table('/Users/sribindusreepada/Desktop/DFCI_CPD/proteomics_viz-master/test/data/BioPlex_293T_Network_10K_Dec_2019.tsv')

def get_hits(og_data, data_hits):
    only_hits = data_hits[["gene", "freq_hit"]]
    only_hits.rename(columns={"gene": "Gene.Symbol"}, inplace=True)
    merged_df = og_data.merge(only_hits, on='Gene.Symbol', how='left')
    merged_df['dot_size'] = merged_df['freq_hit'] * 4

    return merged_df


def get_compounds(data):
    compounds = dict()
    for file_name in data:
        if "CPD." in file_name:
            list_split = file_name.split("_")
            compound = list_split[7]
            compounds[file_name] = compound
        else:
            pass
    return compounds


compounds = get_compounds(data)


def one_column(df):
    all_cols = []
    for column in df.columns[:]:
        col_list = df[column].tolist()
        all_cols.append(col_list)

    flat_list = []
    for elem in all_cols:
        flat_list.extend(elem)

    return flat_list


# adding -logp column to datasets
for key, value in data.items():
    if "P.Value" in value.columns.tolist():
        log_transform = [math.log(i) * -1 for i in value["P.Value"]]
        value["-logp"] = log_transform
        data[key] = value
    else:
        pass

# adding image url column
for key, value in data.items():
    str_drop = str(key)
    if "ip-cpd" in str_drop:
        list_split = str_drop.split("_")
        compound = list_split[7]
        cpd_list = compound.split(".")
        cmp = cpd_list[-1]
        cmp = 'CPD' + cmp
        cmp_smile_url = get_smile_url(cmp)
        images = [cmp_smile_url for i in range(len(value.index))]
        value["images"] = images
        # print(key,compound,cmp_smile_url)
        data[key] = value
    else:
        pass

app = dash.Dash(external_stylesheets=[dbc.themes.FLATLY])

sidebar = html.Div([
    dbc.Row(
        [
            html.H5('Settings',
                    style={'margin-top': '12px', 'margin-bottom': '12px'})
        ],
        style={"height": "7vh", 'textAlign': 'center'},
        className='bg-primary text-white font-italic'
    ),

    dbc.Row([html.Div([
        html.P('Select a CSV file:', style={'margin-top': '8px', 'margin-bottom': '4px'},
               className='font-weight-bold'),

        dcc.Dropdown(
            id='dropdown',
            options=[{'label': i, 'value': i} for i in list(data.keys()) if "ip-cpd" in i],
            value=list(data.keys())[0]
        ),

        html.P("-logP threshold:", style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),
        dcc.Input(
            id="logP_slider", type="number", placeholder="input logP threshold", value=1,
            min=0, max=20
        ),
        # dcc.Slider(id='logP_slider', min=0, max=20, step=1, value=1),
        html.P("logFC thresholds:", style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),
        dcc.Input(
            id="logFC_neg_slider", type="text", placeholder="negative value only", value=-1,
            min=-100, max=0
        ),
        dcc.Input(
            id="logFC_pos_slider", type="text", placeholder="positive value only", value=1,
            min=0.01, max=100
        ),

        html.P('Select a txt file:', style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),

        dcc.Dropdown(
            id='dropdown2',
            options=[{'label': i, 'value': i} for i in list(txt_data.keys())],
            value=list(txt_data.keys())[0]
        ),

        html.P("Compound: ",
               style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),

        dcc.Dropdown(
            id='compound',
            options=[{'label': i, 'value': i} for i in list(compounds.values())]
        ),

        html.Hr(),

        html.P('Compound Structure: ', style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),

        html.Img(id="image", src='',
                 style={'width': '10vw', 'height': '20vh', 'display': 'inline-block', 'vertical-align': 'top'}),

        html.P('Box plot: ', style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),

        dcc.Graph(
                id='box-plot',
                style={'width': '30vw', 'height': '45vh', 'vertical-align': 'top'}
        ),

    ]),
    ], style={'height': '100vh', 'margin': '8px'})
])

content = html.Div([
    dbc.Row(
        [
            html.H3('Proteomics Data Visualization',
                    style={'margin-top': '4px', 'margin-bottom': '3px'})
        ],
        style={"height": "7vh", 'textAlign': 'center'},
        className='bg-primary text-white font-italic'
    ),
    html.Div([
    dcc.Graph(
        id='volcano_plot',
        style={'width': '65vw', 'height': '75vh', 'vertical-align': 'top'}
    ), ], style={'display': 'flex', 'justify-content': 'center', 'align-items': 'center'}),

    html.H3(children='Hover data', style={'font-size': '30px', 'textAlign': 'center'}),

    html.Div(
        id='hover-data'
    ),

    html.Div(
        id='selected-data'
    )
])

app.layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(sidebar, width=4, className='bg-light'),
                dbc.Col(content, width=8)
            ]
        ),
    ],
    fluid=True
)


@app.callback(
    Output('compound', 'value'),
    Input('dropdown', 'value'),
)
def set_compound_options(dropdown):
    str_drop = str(dropdown)
    list_split = str_drop.split("_")
    compound = list_split[7]
    cpd_list = compound.split(".")
    num = cpd_list[-1]
    return num

@app.callback(
    Output('volcano_plot', 'figure'),
    Input('dropdown', 'value'),
    Input('logFC_neg_slider', 'value'),
    Input('logFC_pos_slider', 'value'),
    Input('logP_slider', 'value'),
    Input('volcano_plot', 'selectedData')
)

def update_graph(dropdown, logFC_neg_slider, logFC_pos_slider, logP_slider, selected_data):
    if dropdown is not None:
        df = data[dropdown]
        df2 = get_hits(df, hits)
        df2['freq_hit'] = df2['freq_hit'].replace(np.nan, 0.5)
        df2['dot_size'] = df2['dot_size'].replace(np.nan, 4)

        color = [0 if v["-logp"] < float(logP_slider) or (v["logFC"] > float(logFC_neg_slider) and v["logFC"] < float(logFC_pos_slider)) else 1 if v["logFC"] > float(logFC_pos_slider) and v["-logp"] > float(logP_slider) else -1 for index, v in df[["logFC", "-logp"]].iterrows()]
        colorscale = [[0, 'blue'], [0.5, 'gray'], [1.0, 'red']]
        opacity = [0.2 if v["-logp"] < float(logP_slider) or (v["logFC"] > float(logFC_neg_slider) and v["logFC"] < float(logFC_pos_slider)) else 1 for index, v in df2[["logFC", "-logp"]].iterrows()]
        fig = px.scatter(df2, x="-logp", y="logFC", size="dot_size",
                         hover_data=["freq_hit", "Gene.Symbol",
                                     "Accession"], custom_data=["images"])
        fig.add_vline(x=logP_slider, line_dash="dash")
        fig.add_hline(y=logFC_pos_slider, line_dash="dash")
        fig.add_hline(y=logFC_neg_slider, line_dash="dash")
        fig.update_traces(marker=dict(color=color, opacity=opacity, colorscale=colorscale))
        fig.update_traces(selector=dict(mode='markers+text'))
        fig.update_layout(
            uirevision='same'
        )
        fig.update_layout(clickmode='event+select')
        fig.update_layout(title='Volcano Plot', xaxis_title="-logp", yaxis_title="logFC",
                          hoverlabel=dict(bgcolor="white", font_size=12, font_family="HelveticaNeue"))
        fig.update_layout(margin=dict(l=20, r=20))


        if selected_data:
            selected_point = selected_data['points'][0]
            selected_uniprot = selected_point["customdata"][-1]

            connected_uniprots = network_df[network_df['UniprotA'] == selected_uniprot]['UniprotB'].tolist()
            connected_uniprots.append(selected_uniprot)

            opacity = [1.0 if acc in connected_uniprots else 0.1 for acc in df2['Accession']]
            color2 = ['green' if acc in connected_uniprots else color for acc, color in
                     zip(df2['Accession'], fig['data'][0]['marker']['color'])]

            fig.update_traces(
                marker=dict(color=color2, opacity=opacity)
            )

        return fig


@app.callback(
    Output('image', 'src'),
    Input('volcano_plot', 'hoverData')
)
def open_url(hoverData):
    if hoverData:
        # print(hoverData)
        return hoverData["points"][0]["customdata"][0]
    else:
        raise PreventUpdate


@app.callback(
    Output('box-plot', 'figure'),
    Input('dropdown2', 'value'),
    Input('volcano_plot', 'selectedData'),
    Input('compound', 'value')
)
def get_data(dropdown2, selectedData, compound):
    if selectedData is None:
        fig = px.box()
        return fig
    else:
        df2 = txt_data[dropdown2]
        accession = selectedData["points"][0]["customdata"][-1]
        accession_df = df2.loc[df2['Master Protein Accessions'] == accession]
        dmso = accession_df[["Abundance F1 Sample DMSO", "Abundance F2 Sample DMSO",
                             "Abundance F3 Sample DMSO", "Abundance F4 Sample DMSO"]]
        colNames = accession_df.columns[accession_df.columns.str.contains(pat=compound)]
        colNames = colNames.tolist()
        new_list = list()
        for item in colNames:
            if "Abundance " in item:
                new_list.append(item)
        cpd = accession_df[new_list]
        all = pd.DataFrame()
        all["dmso"] = one_column(dmso)
        all["cpd"] = one_column(cpd)

        # Melt the data
        melted_data = all.melt(value_vars=['dmso', 'cpd'], var_name='column', value_name='value')

        # Replace 'dmso' with 1 and 'cpd' with 2
        melted_data['column'] = melted_data['column'].replace({'dmso': 1, 'cpd': 2})

        # Create the box plot and scatter points
        box_trace = go.Box(
            x=melted_data['column'],
            y=melted_data['value'],
            name='Box Plot',
            showlegend = False
        )

        scatter_trace = go.Scatter(
            x=melted_data['column'],
            y=melted_data['value'],
            mode='markers',
            marker=dict(size=5),
            name='Scatter Points',
            showlegend=False
        )

        data = [box_trace, scatter_trace]

        # Create the figure
        layout = go.Layout(
            title=f'Uniprot {accession}',
            xaxis=dict(title="Compound", tickvals=[1, 2], ticktext=["DMSO", "CPD"]),
            yaxis=dict(title="Abundance"),
            hoverlabel=dict(bgcolor="white", font_size=12, font_family="HelveticaNeue"),
            boxmode='overlay'  # Set boxmode to overlay
        )

        fig = go.Figure(data=data, layout=layout)

        return fig


@app.callback(
    Output('hover-data', 'children'),
    [Input('volcano_plot', 'hoverData')]
)
def display_hover_data(hoverData):
    return html.Pre(json.dumps(hoverData, indent=2))

if __name__ == '__main__':
    app.run_server(debug=True, port=50030)
