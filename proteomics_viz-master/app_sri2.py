import dash
import math
from dash import dcc, html, Input, Output, State, dash_table
import pandas as pd
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

import plotly.express as px
import json
import numpy as np
import plotly.graph_objects as go
from data_loader import load_data, load_txt, get_compounds, get_experiment, get_txtfile_from_data, get_smile_url
from google_auth import GoogleAuth

# external stylesheet
external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

# Load the data
folder = 'data'
data = load_data(folder)
txt_data = load_txt(folder)
hits = pd.read_csv('pde3a_inh_ip_hits_p0001.csv')
hits = hits[["gene", "freq_hit"]]
network_df = pd.read_table('BioPlex_293T_Network_10K_Dec_2019.tsv')
compounds = get_compounds(data)
experiments = get_experiment(txt_data)
MAX_FILENAME_LENGTH = 100


# Iterate through each row in the DataFrame and create gene column
uni_gene_dict = {}
for filename, df in data.items():
    # Iterate through rows in the data frame
    if 'Accession' in df.columns.tolist():
        temp_dict = df.set_index('Accession').to_dict()['Gene.Symbol']
        uni_gene_dict.update(temp_dict)
    # for index, row in df.iterrows():
    #     if 'Accession' in row.
    #         accession = row['Accession']
    #         gene_symbol = row['Gene.Symbol']

    #         # Check if 'Accession' is not already in the mapping to avoid duplicates
    #         if accession not in uni_gene_dict:
    #             uni_gene_dict[accession] = gene_symbol
    else:
        # get first two columns of dataframe
        temp_dict = df.iloc[:, :2]
        temp_dict.columns = ['Accession', 'Gene.Symbol']
        temp_dict = temp_dict.set_index('Accession').to_dict()['Gene.Symbol']
        uni_gene_dict.update(temp_dict)
        

# Get compound  column
for file_name, data_frame in data.items():
    cpd_index = file_name.find('CPD')

    if cpd_index != -1:
        seven_digits = file_name[cpd_index + 4: cpd_index + 11]
        compound_code = f'CPD-{seven_digits}'
        cpd_part = file_name[cpd_index:]
        cpd_part = cpd_part[:-4]

    data_frame['CPD'] = compound_code
    data_frame['Exp'] = cpd_part

def get_hits(og_data, data_hits):
    # only_hits = data_hits[["gene", "freq_hit"]]
    data_hits.rename(columns={"gene": "Gene.Symbol"}, inplace=True)
    # check if gene_symbol is in the columns
    if 'Gene.Symbol' not in og_data.columns.tolist():
        columns = og_data.columns.tolist()
        new_columns = ['Accession', 'Gene.Symbol'] + columns[2:]
        og_data.columns = new_columns
    merged_df = og_data.merge(data_hits, on='Gene.Symbol', how='left')
    max_hit = merged_df['freq_hit'].max()
    merged_df['dot_size'] = merged_df['freq_hit'] * 2
    merged_df['dot_size'] = merged_df['dot_size'].apply(lambda x: 40 if x >= 40 else (10 if x <= 10 else x))
    return merged_df


app = dash.Dash(external_stylesheets=[dbc.themes.FLATLY])
auth = GoogleAuth(app)
sidebar_top = html.Div([
    dbc.Row(
        [
            html.H5('Settings',
                    style={'margin-top': '12px', 'margin-bottom': '12px'})
        ],
        style={"height": "5vh", 'textAlign': 'center'},
        className='bg-primary text-white font-italic'
    )
])
searchbar = html.Div([
    dbc.Row([html.Div([
        html.P('Select a Gene/Protein of Interest', style={'margin-top': '8px', 'margin-bottom': '4px'},
               className='font-weight-bold'),
        dcc.Dropdown(
            id='genes',
            options=[],
            value=None,
            placeholder='Select or Search a gene',
            search_value='',
            clearable=True,
            style={"maxWidth": "100%"}
                ),
        html.P('Select a CSV file:', style={'margin-bottom': '4px'},
               className='font-weight-bold'),
        dcc.Dropdown(
            id='datafile',
            options=[],
            value=None,
            placeholder='Select or Search a CSV file',
            search_value='',
            clearable=True,
            style={"maxWidth": "100%"}
        ),
        html.Br()
    ])])
])

sidebar = html.Div([
    dbc.Row([html.Div([

        html.P("-logP threshold: default 3, which is (10^-3)", style={'margin-top': '8px', 'margin-bottom': '4px'},
               className='font-weight-bold'),
        dcc.Input(
            id="logP_slider", type="number", placeholder="input -logP threshold"
        ),
        # dcc.Slider(id='logP_slider', min=0, max=20, step=1, value=1),
        html.P("logFC thresholds:", style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),
        dcc.RangeSlider(
            id="logFC_slider", step=0.1, dots=False, allowCross=False,
            tooltip={"placement": "bottom", "always_visible": True},
        ),
        html.P('Select a txt file:', style={'margin-top': '16px', 'margin-bottom': '4px'},
               className='font-weight-bold'),
        dcc.Dropdown(
            id='txtfile',
            options=[{'label': i, 'value': i} for i in list(txt_data.keys())],
            value=list(txt_data.keys())[0],
            style={"maxWidth": "100%"}
        ),
        html.Hr(),
        html.Img(id="cmpd_image", src='',
                 style={'width': '50%', 'height': 'auto', 'display': 'inline-block', 'vertical-align': 'top'}),
        html.Div(id='compound',
                 style={'display': 'inline-block', 'vertical-align': 'center', 'horizontal-align': 'center',
                        'margin-left': '10px', 'font-size': '20px'}),
        html.P('Box plot: ', style={'margin-top': '8px', 'margin-bottom': '4px'},
               className='font-weight-bold'),

        dcc.Graph(
            id='box-plot',
            style={'width': '100%', 'height': 'auto', 'vertical-align': 'top', 'bgcolor': 'white'}
        ),

    ]),
    ], style={'height': '100%', 'margin': '8px'})
])

content_top = html.Div([
    dbc.Row(
        [
            html.H3('Proteomics Data Visualization',
                    style={'margin-top': '4px', 'margin-bottom': '3px'})
        ],
        style={"height": "5vh", 'textAlign': 'center'},
        className='bg-primary text-white font-italic'
    )
])
content = html.Div([
    html.Div([
        dcc.Tabs(id='tabs', value='compound_tab', children=[
            dcc.Tab(label='Compound', value='compound_tab', children=[
                dcc.Graph(
                    id='volcano_plot',
                    style={'width': '100%', 'height': '60vh', 'vertical-align': 'top', 'margin-top': '8px', 'bgcolor': 'white'}
                ),
            ]),
            dcc.Tab(label='Gene', value='gene_tab', children=[
                dcc.Graph(
                    id='gene_plot',
                    style={'width': '100%', 'height': '60vh', 'vertical-align': 'top', 'margin-top': '8px', 'bgcolor': 'white'}
                ),
            ]),
        ]),
    ]),
    html.Div(
        id='hover-data'
    ),
    html.Div(
        id='selected-data'
    ),
    html.Br(),
    html.Div([
        html.H5('BioPlex Network Data for Selected Point', style={'textAlign': 'center'}),
        dash_table.DataTable(
            id='data-table',
            columns=[],
            data=[],
            page_size=10,
            sort_action="native",
            filter_action="native",
            style_table={'overflowX': 'auto'},
            style_cell={'textAlign': 'left'},
            style_data={
                'color': 'black',
                'backgroundColor': 'white'
            },
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(220, 220, 220)'
                }
            ]
        )
    ], style={'width': '100%'})
])


app.layout = dbc.Container(
    [
        dbc.Row([
            dbc.Col(sidebar_top, width=4, className='bg-light'),
            dbc.Col(content_top, width=8)
        ]),
        dbc.Row([
            dbc.Col(searchbar, width=12, className='bg-light')
        ]),
        dbc.Row([
            dbc.Col(sidebar, width=4, className='bg-light'),
            dbc.Col(content, width=8)
        ])
    ],
    fluid=True
)

@app.callback(
    Output('genes', 'options'),
    Input('genes', 'search_value')
)
def update_file_options(search_value):
    genes_list = list(uni_gene_dict.values())
    return genes_list

@app.callback(
    Output('datafile', 'options'),
    Input('datafile', 'search_value')
)
def update_file_options(search_value):
    all_files = list(data.keys())
    if not search_value:
        return [
            {'label': (f if len(f) <= MAX_FILENAME_LENGTH else f"{f[:MAX_FILENAME_LENGTH]}..."), 'value': f, 'title': f}
            for f in all_files]
    filtered_files = [f for f in all_files if search_value.lower() in f.lower()]
    if not filtered_files:
        return [{'label': 'No files found', 'value': None, 'title': None}]

    options = []
    for fn in filtered_files:
        label = fn
        if len(fn) > MAX_FILENAME_LENGTH:
            label = f"{fn[:MAX_FILENAME_LENGTH]}..."
        options.append({'label': label, 'value': fn, 'title': fn})
    return options

@app.callback(
    Output('volcano_plot', 'figure'),
    Input('datafile', 'value'),
    Input('logFC_slider', 'value'),
    Input('logP_slider', 'value'),
    Input('volcano_plot', 'selectedData')
)
def update_graph(datafile, logFC_slider_value, logP_slider, selected_data):
    if datafile is None:
        return go.Figure()
    logP_slider = logP_slider if logP_slider is not None else 3
    df_data = data[datafile]
    df_data_hits = get_hits(df_data, hits)
    df_data_hits['freq_hit'] = df_data_hits['freq_hit'].replace(np.nan, 0)
    df_data_hits['dot_size'] = df_data_hits['dot_size'].replace(np.nan, 5)
    logFC_neg_slider, logFC_pos_slider = logFC_slider_value
    conditions = [
        (df_data_hits['logFC'] >= logFC_pos_slider) & (df_data_hits['-logp'] >= logP_slider),
        (df_data_hits['logFC'] <= logFC_neg_slider) & (df_data_hits['-logp'] >= logP_slider),
    ]

    df_data_hits['color'] = np.select(conditions, ['red', 'blue'], 'gray')
    df_data_hits['opacity'] = np.select(conditions, [1, 1], 0.25)
    # print(df_data_hits['color'])
    # print(df_data_hits['opacity'])
    # df_data_hits['color'],df_data_hits['opacity'] = zip(*df_data_hits.apply(lambda row: assign_color_opacity(row,logFC_pos_slider,logFC_neg_slider,logP_slider), axis=1))
    # opacity = [0.25 if v["-logp"] < float(logP_slider) or (v["logFC"] > float(logFC_neg_slider) and v["logFC"] < float(logFC_pos_slider)) else 1 for index, v in df_data_hits[["logFC", "-logp"]].iterrows()]
    fig = px.scatter(df_data_hits, x="-logp", y="logFC", size="dot_size",
                     hover_data=["freq_hit", "Gene.Symbol",
                                 "Accession"])
    fig.add_vline(x=logP_slider, line_dash="dash", line_width=0.5)
    fig.add_hline(y=logFC_pos_slider, line_dash="dash", line_width=0.5)
    fig.add_hline(y=logFC_neg_slider, line_dash="dash", line_width=0.5)
    fig.update_layout(plot_bgcolor='white')
    fig.update_traces(marker=dict(color=df_data_hits['color'], opacity=df_data_hits['opacity']))
    fig.update_traces(selector=dict(mode='markers'))
    # get dataframe meeting conditions
    df_annotation = df_data_hits[(df_data_hits['-logp'] >= logP_slider) & (
                (df_data_hits['logFC'] >= logFC_pos_slider) | (df_data_hits['logFC'] <= logFC_neg_slider))]
    # add annotations for the selected points
    annotations = [{
        'x': row['-logp'],
        'y': row['logFC'],
        'xref': 'x',
        'yref': 'y',
        'text': row['Gene.Symbol'],
        'showarrow': True,
        'font': dict(color='black', size=8)
    } for _, row in df_annotation.iterrows()]

    fig.update_layout(
        uirevision='same'
    )
    fig.update_layout(clickmode='event+select')
    fig.update_layout(title='Volcano Plot', title_x=0.5, title_font=dict(size=24),
                      xaxis=dict(title="-logp", titlefont=dict(size=24), tickfont=dict(size=14)),
                      yaxis=dict(title="logFC", titlefont=dict(size=24), tickfont=dict(size=14)),
                      hoverlabel=dict(bgcolor="white", font_size=16, font_family="HelveticaNeue"))
    fig.update_layout(margin=dict(l=20, r=20))
    fig.update_layout(annotations=annotations)
    if selected_data:
        selected_point = selected_data['points'][0]
        selected_uniprot = selected_point["customdata"][2]
        connected_uniprots_A = network_df[network_df['UniprotA'] == selected_uniprot]['UniprotB'].tolist()
        connected_uniprots_B = network_df[network_df['UniprotB'] == selected_uniprot]['UniprotA'].tolist()
        connected_uniprots = connected_uniprots_A + connected_uniprots_B
        connected_uniprots.append(selected_uniprot)

        # opacity = [1.0 if acc in connected_uniprots else 0.1 for acc in df_data_hits['Accession']]
        # color2 = ['green' if acc in connected_uniprots else color for acc, color in zip(df_data_hits['Accession'], fig['data'][0]['marker']['color'])]

        # Set all points in connected_uniprots to green and 100% opacity
        df_data_hits.loc[df_data_hits['Accession'].isin(connected_uniprots), 'color'] = 'green'
        df_data_hits.loc[df_data_hits['Accession'].isin(connected_uniprots), 'opacity'] = 1.0

        # Set the color and opacity of all other points to their original values
        df_data_hits.loc[~df_data_hits['Accession'].isin(connected_uniprots), 'color'] = 'gray'
        df_data_hits.loc[~df_data_hits['Accession'].isin(connected_uniprots), 'opacity'] = 0.6

        # Update all traces with the modified color and opacity
        fig.update_traces(
            marker=dict(color=df_data_hits['color'], opacity=df_data_hits['opacity']),
            selector=dict(mode='markers')
        )

    else:
        df_data_hits['color'] = np.select(conditions, ['red', 'blue'], 'gray')
        df_data_hits['opacity'] = np.select(conditions, [1, 1], 0.3)
        fig.update_traces(
            marker=dict(color=df_data_hits['color'], opacity=df_data_hits['opacity'])
        )
        # fig.update_traces(
        #     marker=dict(color=df_data_hits['color'], opacity=df_data_hits['opacity'])
        # )
        '''
        for idx in connected_uniprots:
            fig.data[0].marker.color[idx] = 'green'
        '''

    return fig

@app.callback(
    Output('gene_plot', 'figure'),
    Input('genes', 'value'),
    Input('logFC_slider', 'value'),
    Input('logP_slider', 'value')
)

def update_gene_graph(gene, logFC_slider_value, logP_slider):
    if gene is None:
        return go.Figure()
    logP_slider = logP_slider if logP_slider is not None else 3
    all_files = list(data.keys())
    filtered_gene = [f for f in all_files if gene in data[f]["Gene.Symbol"].values]
    filtered_data = {filename: data[filename] for filename in filtered_gene if filename in data}

    filtered_dfs = []
    for filename, df in filtered_data.items():
        gene_rows = df[df['Gene.Symbol'] == gene]
        filtered_dfs.append(gene_rows)

    df_data = pd.concat(filtered_dfs, ignore_index=True)

    logFC_neg_slider, logFC_pos_slider = logFC_slider_value
    conditions = [
        (df_data['logFC'] >= logFC_pos_slider) & (df_data['-logp'] >= logP_slider),
        (df_data['logFC'] <= logFC_neg_slider) & (df_data['-logp'] >= logP_slider),
    ]

    df_data['color'] = np.select(conditions, ['red', 'blue'], 'gray')
    df_data['opacity'] = np.select(conditions, [1, 1], 0.25)
    # print(df_data_hits['color'])
    # print(df_data_hits['opacity'])
    # df_data_hits['color'],df_data_hits['opacity'] = zip(*df_data_hits.apply(lambda row: assign_color_opacity(row,logFC_pos_slider,logFC_neg_slider,logP_slider), axis=1))
    # opacity = [0.25 if v["-logp"] < float(logP_slider) or (v["logFC"] > float(logFC_neg_slider) and v["logFC"] < float(logFC_pos_slider)) else 1 for index, v in df_data_hits[["logFC", "-logp"]].iterrows()]
    fig = px.scatter(df_data, x="-logp", y="logFC",
                     hover_data=["CPD", "Exp", "Gene.Symbol",
                                 "Accession"])
    fig.add_vline(x=logP_slider, line_dash="dash", line_width=0.5)
    fig.add_hline(y=logFC_pos_slider, line_dash="dash", line_width=0.5)
    fig.add_hline(y=logFC_neg_slider, line_dash="dash", line_width=0.5)
    fig.update_layout(plot_bgcolor='white')
    fig.update_traces(marker=dict(color=df_data['color'], opacity=df_data['opacity']))
    fig.update_traces(selector=dict(mode='markers'))
    # get dataframe meeting conditions
    df_annotation = df_data[(df_data['-logp'] >= logP_slider) & (
                (df_data['logFC'] >= logFC_pos_slider) | (df_data['logFC'] <= logFC_neg_slider))]
    # add annotations for the selected points
    annotations = [{
        'x': row['-logp'],
        'y': row['logFC'],
        'xref': 'x',
        'yref': 'y',
        'text': row['CPD'],
        'showarrow': True,
        'font': dict(color='black', size=8)
    } for _, row in df_annotation.iterrows()]

    fig.update_layout(
        uirevision='same'
    )
    fig.update_layout(clickmode='event+select')
    fig.update_layout(title='Volcano Plot', title_x=0.5, title_font=dict(size=24),
                      xaxis=dict(title="-logp", titlefont=dict(size=24), tickfont=dict(size=14)),
                      yaxis=dict(title="logFC", titlefont=dict(size=24), tickfont=dict(size=14)),
                      hoverlabel=dict(bgcolor="white", font_size=16, font_family="HelveticaNeue"))
    fig.update_layout(margin=dict(l=20, r=20))
    fig.update_layout(annotations=annotations)

    return fig

@app.callback(
    [Output('compound', 'children'),
     Output('logFC_slider', 'min'),
     Output('logFC_slider', 'max'),
     Output('logFC_slider', 'value'),
     Output('logFC_slider', 'marks'),
     Output('txtfile', 'options'),
     Output('txtfile', 'value'),
     Output('cmpd_image', 'src')],
    [Input('datafile', 'value'),
     Input('tabs', 'value'),
     Input('gene_plot', 'selectedData')]
)
def set_compound_options(datafile, selected_tab, selected_data):
    if datafile:
        data_value = datafile
        cmpd = data_value.split('output_')[1].split('_vs')[0]
        compound = cmpd.replace('.', '-')

        # compound = compounds[data_value]
        df = data[data_value]
        data_min_value = df['logFC'].min()
        data_max_value = df['logFC'].max()
        range_min = math.floor(data_min_value)
        range_max = math.ceil(data_max_value)
        data_median = df['logFC'].median()
        logfc_slider_values = [data_median - 1, data_median + 1]
        logfc_marks = {i: str(i) for i in range(range_min, range_max + 1)}
        compound_url = get_smile_url(compounds[data_value])
        # print(data_value,cmpd,compound,compound_url)
        txtfile = get_txtfile_from_data(data_value, txt_data)
        # print(txtfile)

        if selected_tab == 'gene_tab':
            if selected_data:
                selected_point = selected_data['points'][0]
                compound = selected_point['customdata'][0]
                compound_url = get_smile_url(compound)

    else:
        compound = ""
        range_min = -2
        range_max = 2
        logfc_slider_values = [-1.5, 1.5]
        logfc_marks = {i: str(i) for i in range(range_min, range_max + 1)}
        txtfile = list(txt_data.keys())
        compound_url = ""
    return compound, range_min, range_max, logfc_slider_values, logfc_marks, [{'label': i, 'value': i} for i in
                                                                              txtfile], txtfile[0], compound_url

@app.callback(
    Output('box-plot', 'figure'),
    Input('txtfile', 'value'),
    Input('volcano_plot', 'selectedData'),
    Input('compound', 'children')
)
def get_data(txtfile, selectedData, children):
    if selectedData is None:
        fig = px.box()
        return fig
    else:
        df_txt = txt_data[txtfile]
        accession = selectedData["points"][0]
        accession = accession['customdata'][2]
        gene = selectedData['points'][0]['customdata'][1]
        accession_df = df_txt.loc[df_txt['Master Protein Accessions'] == accession]
        accession_df_columns = accession_df.columns.tolist()
        # print(accession_df_columns)
        # get the columns that contains dmso or DMSO
        cmpd = children
        # count '-' in cmpd
        # print(accession_df_columns,cmpd)
        tag = cmpd.split('-')[-1] if cmpd.count('-') >= 2 else None
        dmso_tag = f'dmso-{tag}' if tag else 'dmso'
        dmso_tag = dmso_tag.lower()
        dmso_cols = [col for col in accession_df_columns if
                     dmso_tag in col.lower() and 'count' not in col.lower() and 'abundance' in col.lower()]
        df_dmso = accession_df[dmso_cols]
        cmpd_columns = [col for col in accession_df_columns if
                        cmpd in col and 'count' not in col.lower() and 'abundance' in col.lower()]
        df_cpd = accession_df[cmpd_columns]
        dmso_values = df_dmso.values.flatten().tolist()
        cpd_values = df_cpd.values.flatten().tolist()
        # dmso_values = [x for x in dmso_values if str(x) != 'nan']
        # cpd_values = [x for x in cpd_values if str(x) != 'nan']
        box_data = [(dmso_tag.upper(), x) for x in dmso_values] + [(cmpd, x) for x in cpd_values]
        df_box = pd.DataFrame(box_data, columns=['compound', 'value'])
        # print(df_dmso.columns.tolist(),df_cpd.columns.tolist())
        # print(df_box)
        fig = px.strip(df_box, x='compound', y='value', color='compound', hover_data=['compound', 'value'],
                       stripmode='overlay')
        fig.add_trace(
            go.Box(y=df_box.query("compound == @dmso_tag.upper()")['value'], name=dmso_tag.upper(), boxpoints=False))
        fig.add_trace(go.Box(y=df_box.query("compound == @cmpd")['value'], name=cmpd, boxpoints=False))
        fig.update_layout(title=f'Uniprot:{accession}/Gene:{gene}', title_x=0.5, xaxis_title='Compound',
                          yaxis_title='Abundance', plot_bgcolor='white',
                          hoverlabel=dict(bgcolor="white", font_size=12, font_family="HelveticaNeue"), showlegend=False)
        return fig


@app.callback(
    Output('hover-data', 'children'),
    Input('volcano_plot', 'hoverData')
)
def display_hover_data(hoverData):
    if hoverData is None:
        hover = ""
    else:
        hover_point = hoverData["points"][0]
        x = hover_point["x"]
        y = hover_point["y"]
        customdata = hover_point["customdata"]
        # print(customdata)
        hit, gene, uniprot = customdata[:3]
        hover = f"Hover Data: Gene: {gene}, Uniprot: {uniprot}, -logp: {x}, logFC: {y}, hits: {hit}"
    return hover


@app.callback(
    [Output('selected-data', 'children'),
     Output('data-table', 'columns'),
     Output('data-table', 'data')],
    Input('gene_plot', 'selectedData')
)
def display_select_data(selectData):
    # print(selectData)
    if selectData is None:
        selected = ""
        data_columns = []
        data_table = []
    else:
        select_point = selectData["points"][0]
        x_s = select_point["x"]
        y_s = select_point["y"]
        customdata_s = select_point["customdata"]
        # print(customdata_s)
        hit_s, gene_s, uniprot_s = customdata_s[:3]
        connected_uniprots = network_df[(network_df['UniprotA'] == uniprot_s) | (network_df['UniprotB'] == uniprot_s)]
        selected = f"Selected Data: Gene: {gene_s}, Uniprot: {uniprot_s}, -logp: {x_s}, logFC: {y_s}, hits: {hit_s}"
        if connected_uniprots is not None:
            # keep pW, pNI, pInt with two decimal places in scientific notation
            connected_uniprots['pW'] = connected_uniprots['pW'].apply(lambda x: '%.3E' % x)
            connected_uniprots['pNI'] = connected_uniprots['pNI'].apply(lambda x: '%.3E' % x)
            connected_uniprots['pInt'] = connected_uniprots['pInt'].apply(lambda x: '%.3E' % x)
            data_table = connected_uniprots.to_dict('records')
            data_columns = [{"name": i, "id": i} for i in connected_uniprots.columns]
        else:
            data_table = []
            data_columns = []
    return selected, data_columns, data_table


if __name__ == '__main__':
    app.run_server(host='0.0.0.0', debug=True, port=5000)
