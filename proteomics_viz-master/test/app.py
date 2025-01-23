import dash
import math
from dash import dcc, html, Input, Output, State
import pandas as pd
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import plotly.express as px
import json
import numpy as np
import plotly.graph_objects as go
from data_loader import load_data, load_txt, get_compounds,get_experiment, get_txtfile_from_data,get_smile_url
from google_auth import GoogleAuth
# external stylesheet
external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

# Load the data
folder = 'data'
data = load_data(folder)
txt_data = load_txt(folder)
hits = pd.read_csv('pde3a_inh_ip_hits_p0001.csv')
network_df = pd.read_table('BioPlex_293T_Network_10K_Dec_2019.tsv')
compounds = get_compounds(data)
experiments = get_experiment(txt_data)

def get_hits(og_data, data_hits):
    only_hits = data_hits[["gene", "freq_hit"]]
    only_hits.rename(columns={"gene": "Gene.Symbol"}, inplace=True)
    merged_df = og_data.merge(only_hits, on='Gene.Symbol', how='left')
    merged_df['dot_size'] = merged_df['freq_hit'] * 4
    return merged_df


def one_column(df):
    all_cols = []
    for column in df.columns[:]:
        col_list = df[column].tolist()
        all_cols.append(col_list)

    flat_list = []
    for elem in all_cols:
        flat_list.extend(elem)

    return flat_list

app = dash.Dash(external_stylesheets=[dbc.themes.FLATLY])
# auth = GoogleAuth(app)
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
            id='datafile',
            options=[{'label': i, 'value': i} for i in list(data.keys())],
            value=list(data.keys())[0]
        ),
        html.P("-logP threshold:", style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),
        dcc.Input(
            id="logP_slider", type="number", placeholder="input logP threshold",value=1
        ),
        # dcc.Slider(id='logP_slider', min=0, max=20, step=1, value=1),
        html.P("logFC thresholds:", style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),
        dcc.RangeSlider(
            id="logFC_slider", step = 0.01, dots = False, allowCross=False, tooltip={"placement": "bottom","always_visible": True},
        ),
        html.P('Select a txt file:', style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),
        dcc.Dropdown(
            id='txtfile',
            options=[{'label': i, 'value': i} for i in list(txt_data.keys())],
            value=list(txt_data.keys())[0]
        ),
        html.Hr(),
        html.P('Compound Structure: ', style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),
        html.Img(id="cmpd_image", src='',
                 style={'width': '10vw', 'height': '20vh', 'display': 'inline-block', 'vertical-align': 'top'}),
        html.Div(id='compound'),
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
    [Output('compound', 'children'),
     Output('logFC_slider','min'),
     Output('logFC_slider','max'),
     Output('logFC_slider','value'),
     Output('logFC_slider','marks'),
     Output('txtfile', 'options'),
     Output('txtfile','value'),
     Output('cmpd_image','src')],
    Input('datafile', 'value'),
)
def set_compound_options(datafile):
    data_value = datafile if datafile else list(data.keys())[0]
    # data_value = str(datafile)
    compound = compounds[data_value]
    df = data[data_value]
    data_min_value =df['logFC'].min()
    data_max_value = df['logFC'].max()
    range_min = math.floor(data_min_value)
    range_max = math.ceil(data_max_value)
    data_median = df['logFC'].median()
    logfc_slider_values = [data_median-1,data_median+1]
    logfc_marks = {i: str(i) for i in range(range_min,range_max+1)}
    compound_url = get_smile_url(compound)
    txtfile = get_txtfile_from_data(datafile,txt_data)
    return compound,range_min,range_max,logfc_slider_values,logfc_marks,[{'label': i, 'value': i} for i in txtfile],txtfile[0],compound_url


@app.callback(
    Output('volcano_plot', 'figure'),
    Input('datafile', 'value'),
    Input('logFC_slider', 'value'),
    Input('logP_slider', 'value'),
    Input('volcano_plot', 'selectedData')
)

def update_graph(datafile, logFC_slider_value, logP_slider, selected_data):
    if datafile is not None:
        df = data[datafile]
        df2 = get_hits(df, hits)
        df2['freq_hit'] = df2['freq_hit'].replace(np.nan, 0.5)
        df2['dot_size'] = df2['dot_size'].replace(np.nan, 4)
        logFC_neg_slider, logFC_pos_slider = logFC_slider_value
        color = [0 if v["-logp"] < float(logP_slider) or (v["logFC"] > float(logFC_neg_slider) and v["logFC"] < float(logFC_pos_slider)) else 1 if v["logFC"] > float(logFC_pos_slider) and v["-logp"] > float(logP_slider) else -1 for index, v in df[["logFC", "-logp"]].iterrows()]
        colorscale = [[0, 'blue'], [0.5, 'gray'], [1.0, 'red']]
        opacity = [0.2 if v["-logp"] < float(logP_slider) or (v["logFC"] > float(logFC_neg_slider) and v["logFC"] < float(logFC_pos_slider)) else 1 for index, v in df2[["logFC", "-logp"]].iterrows()]
        fig = px.scatter(df2, x="-logp", y="logFC", size="dot_size",
                         hover_data=["freq_hit", "Gene.Symbol",
                                     "Accession"])
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
            print()
            connected_uniprots = network_df[network_df['UniprotA'] == selected_uniprot]['UniprotB'].tolist()
            connected_uniprots.append(selected_uniprot)

            # opacity = [1.0 if acc in connected_uniprots else 0.1 for acc in df2['Accession']]
            color = ['green' if acc in connected_uniprots else color for acc, color in
                     zip(df2['Accession'], fig['data'][0]['marker']['color'])]

            fig.update_traces(
                marker=dict(color=color, opacity=1.0)
            )

            '''
            for idx in connected_uniprots:
                fig.data[0].marker.color[idx] = 'green'
            '''

        return fig

@app.callback(
    Output('box-plot', 'figure'),
    Input('txtfile', 'value'),
    Input('volcano_plot', 'selectedData'),
    Input('datafile', 'value')
)
def get_data(txtfile, selectedData, datafile):
    if selectedData is None:
        fig = px.box()
        return fig
    else:
        df2 = txt_data[txtfile]
        accession = selectedData["points"][0]
        accession = accession['customdata'][-1]
        accession_df = df2.loc[df2['Master Protein Accessions'] == accession]
        accession_df_columns = accession_df.columns.tolist()
        # print(accession_df_columns)
        # get the columns that contains dmso or DMSO
        dmso_cols = [col for col in accession_df_columns if 'dmso' in col.lower() and 'count' not in col.lower() and 'abundance' in col.lower()]
        # print(dmso_cols)
        df_dmso = accession_df[dmso_cols]
        cmpd = datafile.split('output_')[1].split('_vs')[0]
        cmpd = cmpd.replace('.', '-')
        cmpd_columns = [col for col in accession_df_columns if cmpd in col and 'count' not in col.lower() and 'abundance' in col.lower()]
        # print(cmpd_columns)
        df_cpd = accession_df[cmpd_columns]
        all = pd.DataFrame()
        all["dmso"] = one_column(df_dmso)
        all[cmpd] = one_column(df_cpd)

        # Melt the data
        melted_data = all.melt(value_vars=['dmso', cmpd], var_name='column', value_name='value')

        # Replace 'dmso' with 1 and 'cpd' with 2
        melted_data['column'] = melted_data['column'].replace({'dmso': 1, cmpd: 2})

        # Create the box plot and scatter points
        box_trace = go.Box(
            x=melted_data['column'],
            y=melted_data['value'],
            name='Box Plot',
            showlegend=False
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
            xaxis=dict(title="Compound", tickvals=[1, 2], ticktext=["DMSO", cmpd]),
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
    app.run_server(host='0.0.0.0',debug=True, port=5000)
