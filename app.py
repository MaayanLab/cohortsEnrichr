import os
import json
import pandas as pd
import numpy as np
import dash
import dash_auth
import dash_table as dt
import dash_core_components as dcc
import dash_html_components as html
from react_scatter_board import DashScatterBoard
from scipy.stats import zscore
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate

from dotenv import load_dotenv
load_dotenv()

def yaml_loads(s):
    import io, yaml
    fr = io.StringIO(s)
    return yaml.load(fr)

def yaml_dumps(obj):
    import io, yaml
    fw = io.StringIO()
    yaml.dump(obj, fw)
    return fw.getvalue()

DATA = os.environ.get('DATA', os.path.join(os.path.dirname(__file__), 'data'))
df = pd.read_csv(
    os.path.join(DATA, 'df.tsv'),
    sep='\t'
)
df.index = df.index.astype(str)
df.columns = df.columns.astype(str)
df_umap = pd.read_csv(
    os.path.join(DATA, 'df_umap.tsv'),
    sep='\t',
    index_col='Barcode',
)
df_umap.index = df_umap.index.astype(str)
df_umap.columns = df_umap.columns.astype(str)
df_enrich = pd.read_csv(
    os.path.join(DATA, 'df_enrich.tsv'),
    sep='\t'
)
df_enrich.index = df_enrich.index.astype(str)
df_enrich.columns = df_enrich.columns.astype(str)
df_metadata = pd.read_csv(
    os.path.join(DATA, 'metadata.csv'),
    sep=',',
    index_col='Barcode',
)
df_metadata.index = df_metadata.index.astype(str)
df_metadata.columns = df_metadata.columns.astype(str)
df_cluster_aucs = pd.read_csv(
    os.path.join(DATA, 'cluster_aucs.csv'),
    sep=',',
    index_col=0,
)
df_cluster_aucs.index = df_cluster_aucs.index.astype(str)
df_cluster_aucs.columns = df_cluster_aucs.columns.astype(str)
# df_cluster_aucs.loc[:,:] = zscore(df_cluster_aucs)

meta_cols = [
    'Barcode',
    'Cluster',
    *df_metadata.columns
]

def figure(Barcode=None):
    data = [
        dict(
            Barcode=barcode,
            x=record['UMAP-1'],
            y=record['UMAP-2'],
            label=yaml_dumps(
                dict(
                    Barcode=str(barcode),
                    Cluster=int(record['Cluster']),
                )
            ).replace('\n', '<br>'),
            **record.to_dict(),
        )
        for barcode, record in pd.merge(left=df_umap, left_index=True, right=df_metadata, right_index=True).iterrows()
    ]
    return data


app = dash.Dash(
    __name__,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width"}
    ],
    routes_pathname_prefix=os.environ.get('PREFIX', '')
)
auth = dash_auth.BasicAuth(
    app,
    json.loads(os.environ.get('CREDENTIALS', '{"admin":"admin"}'))
)

app.index_string = '''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css" integrity="sha384-Vkoo8x4CGsO3+Hhxv8T/Q5PaXtkKtu6ug5TOeNV6gBiFeWPGFN9MuhOf23Q9Ifjh" crossorigin="anonymous">
    </head>
    <body>
        <div class="row">
            <div class="col-sm-12 offset-lg-1 col-lg-10">
                <div class="row">
                    <div class="col-sm-12">
                        <h1>cohortsEnrichr</h1>
                        <hr />
                    </div>
                    <div class="col-sm-12">
                        {%app_entry%}
                    </div>
                </div>
            </div>
        </div>
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
'''

app.layout = html.Div(className='row', children=[
    html.Div(
        className='col-sm-8',
        children=[
            DashScatterBoard(
                id='umap',
                data=figure(),
                shapeKey='Cluster',
                colorKey='Cluster',
                labelKeys=['Barcode', 'Cluster'],
                searchKeys=[col for col in df_metadata.columns if col not in ['orig_id']],
                width=800,
                height=500,
                is3d=False,
            ),
        ],
    ),
    html.Div(
        className='col-sm-4',
        children=[
            dt.DataTable(
                id='metadata-table',
                columns=[
                    { 'name': 'attribute', 'id': 'attribute' },
                    { 'name': 'value', 'id': 'value' },
                ],
                sort_action='native',
                sort_mode='multi',
                filter_action='native',
                page_action='native',
                style_header={
                    'backgroundColor': 'rgb(200, 200, 200)',
                    'fontWeight': 'bold'
                },
                style_data={
                    'whiteSpace': 'normal',
                    'height': 'auto',
                },
                style_table={
                    'overflow': 'auto',
                    'height': 450,
                    'paddingRight': 15,
                    'paddingLeft': 15,
                },
                css=[
                    {
                        'selector': '.dash-cell div.dash-cell-value',
                        'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;',
                    },
                ],
                style_data_conditional=[
                    {
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgb(230, 230, 230)'
                    }
                ],
            )
        ]
    ),
    html.Div(
        className='col-sm-12',
        children=[
            html.H2(id='cluster-header'),
        ],
    ),
    html.Div(
        className='col-sm-6',
        children=[
            html.H3('Genetic Enrichment'),
            html.P('We perform cluster vs rest differential expression for each cluster, submit the most significant genesets to Enrichr, and highlight the top enriched terms here. Full Enrichr results link below.'),
            html.Label(id='enrichr-link'),
        ],
    ),
    html.Div(
        className='col-sm-6',
        children=[
            html.H3('Clinical Predictors'),
            html.P('We fit a Logistic Regression using only the attribute in question in an attempt to classify membership in a specific cluster. The AUC of the resulting classifier is recorded, and the relative AUCs reported as Z-Scores.'),
        ],
    ),
    html.Div(
        className='col-sm-6',
        children=[
            dt.DataTable(
                id='data-table',
                columns=[
                    {'name': 'rank', 'id': 'rank'},
                    {'name': 'direction', 'id': 'direction'},
                    {'name': 'term', 'id': 'term'},
                    {'name': 'category', 'id': 'category'},
                    {'name': 'pvalue', 'id': 'pvalue', 'type': 'numeric', 'format': { 'specifier': '.3' } },
                    {'name': 'library', 'id': 'library'},
                ],
                sort_action='native',
                sort_mode='multi',
                filter_action='native',
                filter_query='{pvalue} < 0.05 && {direction} = up',
                page_action='native',
                sort_by=[{ 'column_id': 'pvalue', 'direction': 'asc' }],
                style_as_list_view=True,
                style_header={
                    'backgroundColor': 'rgb(200, 200, 200)',
                    'fontWeight': 'bold'
                },
                style_table={
                    'overflow': 'auto',
                    'width': '100%',
                    'minWidth': '100%',
                    'paddingRight': 15,
                    'paddingLeft': 15,
                },
                style_data={
                    'whiteSpace': 'normal',
                    'height': 'auto',
                },
                css=[
                    {
                        'selector': '.dash-cell div.dash-cell-value',
                        'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;',
                    },
                ],
                style_data_conditional=[
                    {
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgb(230, 230, 230)'
                    }
                ],
                style_cell_conditional=[
                    {
                        'if': { 'column_id': 'term' },
                        'textAlign': 'left',
                        'maxWidth': '20vw',
                    },
                    {
                        'if': { 'column_id': 'category' },
                        'textAlign': 'left',
                    },
                    {
                        'if': { 'column_id': 'direction' },
                        'textAlign': 'center',
                    },
                    {
                        'if': { 'column_id': 'library' },
                        'textAlign': 'left',
                        'maxWidth': '10vw',
                    },
                ],
            ),
        ],
    ),
    html.Div(
        className='col-sm-6',
        children=[
            dt.DataTable(
                id='summary-table',
                columns=[
                    { 'name': 'attribute', 'id': 'attribute' },
                ] + [
                    {'name': col, 'id': col, 'type': 'numeric', 'format': {'specifier': '.3'} }
                    for col in df_cluster_aucs.columns
                ],
                sort_action='native',
                sort_mode='multi',
                filter_action='native',
                page_action='native',
                style_as_list_view=True,
                style_header={
                    'backgroundColor': 'rgb(200, 200, 200)',
                    'fontWeight': 'bold'
                },
                style_table={
                    'overflow': 'auto',
                    'width': '100%',
                    'minWidth': '100%',
                    'paddingRight': 15,
                    'paddingLeft': 15,
                },
                style_data={
                    'whiteSpace': 'normal',
                    'height': 'auto',
                },
                css=[
                    {
                        'selector': '.dash-cell div.dash-cell-value',
                        'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;',
                    },
                ],
                style_data_conditional=[
                    {
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgb(230, 230, 230)'
                    }
                ],
            ),
        ],
    ),
])

lock = False
prevClickData = None
clusterData = None
patientData = {}

@app.callback(
    [
        Output('cluster-header', 'children'),
        Output('enrichr-link', 'children'),
        Output('data-table', 'data'),
        Output('metadata-table', 'data'),
        Output('summary-table', 'data'),
        Output('umap', 'data'),
    ],
    [
        Input('umap', 'clickData'),
        Input('umap', 'hoverData'),
    ]
)
def update_click(clickData, hoverData):
    global lock, prevClickData, clusterData, patientData
    # Initial state
    if not clickData and not hoverData:
        return [
            'Click to cluster to select',
            '',
            [],
            [],
            [],
            figure(),
        ]
    # Get relevant evt
    if prevClickData != clickData: # Click
        lock = not lock
        prevClickData = clickData
        evt = clickData
    else:
        evt = hoverData
    # Get point
    pointData = yaml_loads(evt['label'].replace('<br>', '\n'))
    # Get patient data
    if not lock:
        Barcode = pointData['Barcode']
        patientData['Barcode'] = Barcode
    else:
        Barcode = patientData.get('Barcode')
    # Get cluster
    cluster = pointData['Cluster']
    m = dict(
        **pointData,
        **df_metadata.loc[[pointData['Barcode']]].to_dict('records')[0],
    )
    metadata = [
        {
            'attribute': k,
            'value': m[k],
        }
        for k in meta_cols
    ]
    if not lock and (clusterData is None or clusterData['cluster'] != cluster):
        summary = df_cluster_aucs.reset_index().rename({ 'index': 'attribute' }, axis=1).sort_values(str(cluster), ascending=False).to_dict('records')
        matches = df_enrich[df_enrich['cluster'] == cluster]
        clusterData = {
            'cluster': cluster,
            'summary': summary,
            'matches': matches,
        }
    else:
        cluster = clusterData['cluster']
        summary = clusterData['summary']
        matches = clusterData['matches']
    if matches.size == 0:
        return [
            'Cluster {} ({} samples)'.format(cluster, df_umap[df_umap['Cluster'] == cluster].shape[0]),
            'No data for this cluster',
            [],
            metadata,
            summary,
            figure(Barcode),
        ]
    # Update
    link = matches.iloc[0]['link']
    data = matches.to_dict('records')
    return [
        'Cluster {} ({} samples)'.format(cluster, df_umap[df_umap['Cluster'] == cluster].shape[0]),
        ['Enrichr Link for Cluster ', html.A(link, href=link)],
        data,
        metadata,
        summary,
        figure(Barcode),
    ]

if __name__ == "__main__":
    app.run_server(
        host=os.environ.get('HOST', '0.0.0.0'),
        port=json.loads(os.environ.get('PORT', '8050')),
        debug=json.loads(os.environ.get('DEBUG', 'true')),
    )
