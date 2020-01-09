import os
import json
import pandas as pd
import dash
import dash_auth
import dash_table as dt
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate

from dotenv import load_dotenv
load_dotenv()

DATA = os.environ.get('DATA', os.path.join(os.path.dirname(__file__), 'data'))
df = pd.read_csv(
    os.path.join(DATA, 'df.tsv'),
    sep='\t'
)
df_tsne = pd.read_csv(
    os.path.join(DATA, 'df_tsne.tsv'),
    sep='\t'
)
df_enrich = pd.read_csv(
    os.path.join(DATA, 'df_enrich.tsv'),
    sep='\t'
)

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
        <div class="container">
            <div class="row">
                <div class="col-sm-12">
                    <h1>scEnrichr: An Enrichr interface for 10x Genomics</h1>
                </div>
                <div class="col-sm-12">
                    {%app_entry%}
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

app.layout = html.Div(children=[
    dcc.Graph(
        id='tsne',
        figure={
            'data': [
                dict(
                    x=df_tsne[df_tsne['Cluster'] == cluster]['TSNE-1'],
                    y=df_tsne[df_tsne['Cluster'] == cluster]['TSNE-2'],
                    mode='markers',
                    opacity=1,
                    marker=dict(
                        size=8,
                        line=dict(width=0.5, color='white'),
                    ),
                    text=df_tsne[df_tsne['Cluster'] == cluster]['Cluster'],
                    name='Cluster {} '.format(cluster),
                ) for cluster in sorted(df_tsne['Cluster'].unique())
            ],
            'layout': dict(
                xaxis=dict(
                    title='TSNE-1',
                ),
                yaxis=dict(
                    title='TSNE-2',
                ),
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest'
            ),
        }
    ),
    html.H2(id='cluster-header'),
    html.Label(id='enrichr-link'),
    dt.DataTable(
        id='data-table',
        columns=[
            {'name': 'cluster', 'id': 'cluster'},
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
            'overflowX': 'auto',
            'width': '100%',
            'minWidth': '100%',
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
])

lock = False
prevClickData = None

@app.callback(
    [Output('cluster-header', 'children'), Output('enrichr-link', 'children'), Output('data-table', 'data')],
    [Input('tsne', 'clickData'), Input('tsne', 'hoverData')]
)
def update_click(clickData, hoverData):
    global lock, prevClickData
    # Initial state
    if not clickData and not hoverData:
        return ['Click to cluster to select', '', []]
    # Get relevant evt
    if prevClickData != clickData: # Click
        lock = not lock
        prevClickData = clickData
        evt = clickData
    elif lock: # Hover but locked
        raise PreventUpdate
    else: # Hover and not locked
        evt = hoverData
    # Get cluster
    cluster = evt['points'][0]['text']
    matches = df_enrich[df_enrich['cluster'] == cluster]
    if matches.size == 0:
        return [
            'Cluster {} ({} samples)'.format(cluster, df_tsne[df_tsne['Cluster'] == cluster].shape[0]),
            'No data for this cluster',
            [],
        ]
    # Update
    link = matches.iloc[0]['link']
    data = matches.to_dict('records')
    return [
        'Cluster {} ({} samples)'.format(cluster, df_tsne[df_tsne['Cluster'] == cluster].shape[0]),
        ['Enrichr Link for Cluster ', html.A(link, href=link)],
        data,
    ]

if __name__ == "__main__":
    app.run_server(
        host=os.environ.get('HOST', '0.0.0.0'),
        port=json.loads(os.environ.get('PORT', '8050')),
        debug=json.loads(os.environ.get('DEBUG', 'true')),
    )
