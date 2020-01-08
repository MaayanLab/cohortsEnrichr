import os
import json
import pandas as pd
import dash
import dash_auth
import dash_table as dt
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

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

app = dash.Dash(__name__, meta_tags=[
                {"name": "viewport", "content": "width=device-width"}])
auth = dash_auth.BasicAuth(
    app,
    json.loads(os.environ.get('CREDENTIALS', '{"admin":"admin"}'))
)

app.layout = html.Div(children=[
    html.H1(
        children='10X Genomics Enrichment Dashboard'),

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
                    name='Cluster {}'.format(cluster),
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
    html.Label(
        [html.A(id='link')]
    ),
    dt.DataTable(
        id='table',
        columns=[
            {'name': 'cluster', 'id': 'cluster'},
            {'name': 'rank', 'id': 'rank'},
            {'name': 'direction', 'id': 'direction'},
            {'name': 'term', 'id': 'term'},
            {'name': 'category', 'id': 'category'},
            {'name': 'pvalue', 'id': 'pvalue'},
            {'name': 'library', 'id': 'library'},
        ],
        sort_action='native',
        sort_mode='multi',
        filter_action='native',
        filter_query='{pvalue} < 0.05 && {direction} = up',
        row_selectable='multi',
        column_selectable='multi',
        page_action='native',
        sort_by=[dict(column_id='category', direction='asc'),
                 dict(column_id='rank', direction='asc')],
    ),
])


@app.callback(
    [Output('link', 'children'), Output('link', 'href')],
    [Input('tsne', 'hoverData')]
)
def update_link_on(hoverData):
    if not hoverData:
        return ['', '']
    matches = df_enrich[
        df_enrich['cluster'] == hoverData['points'][0]['text']
    ]
    if matches.size == 0:
        return ['', '']
    link = matches.iloc[0]['link']
    return [link, link]


@app.callback(
    Output('table', 'data'),
    [Input('tsne', 'hoverData')]
)
def update_table_on_hover(hoverData):
    if not hoverData:
        return df_enrich.to_dict('records')
    return df_enrich[df_enrich['cluster'] == hoverData['points'][0]['text']].to_dict('records')


if __name__ == "__main__":
    app.run_server(
        host=os.environ.get('HOST', '0.0.0.0'),
        port=json.loads(os.environ.get('PORT', '8050')),
        debug=json.loads(os.environ.get('DEBUG', 'true')),
    )
