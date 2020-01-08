import sys
import pandas as pd
from collections import OrderedDict

n_genes = 250
top_n_results = 5
# transcription: archs4 coexpression, encode & chEA
# pathways: wikipathways ()
# cell types: human gene atlas, mouse gene atlas, archs4 tissues
useful_libs = OrderedDict([
  ('cell_type', ['Human_Gene_Atlas', 'Mouse_Gene_Atlas', 'ARCHS4_Tissues']),
  ('pathways', ['WikiPathways_2019_Mouse', 'WikiPathways_2019_Human']),
  ('transcription', ['ARCHS4_TFs_Coexp', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X']),
])

# Load data
base_path = sys.argv[1]

df = pd.read_csv(base_path + '/diffexp/graphclust/differential_expression.csv')
df_tsne = pd.read_csv(base_path + '/tsne/2_components/projection.csv')
df_pca = pd.read_csv(base_path + '/pca/10_components/projection.csv')
df_clusters = pd.read_csv(base_path + '/clustering/graphclust/clusters.csv')
df_clusters['Cluster'] = df_clusters['Cluster'].astype(str)

# Util functions
def enrichr_link_from_genes(genes, description='', enrichr_link='https://amp.pharm.mssm.edu/Enrichr'):
  ''' Functional access to Enrichr API
  '''
  import time, requests
  time.sleep(1)
  resp = requests.post(enrichr_link + '/addList', files={
    'list': (None, '\n'.join(genes)),
    'description': (None, description),
  })
  if resp.status_code != 200:
    raise Exception('Enrichr failed with status {}: {}'.format(
      resp.status_code,
      resp.text,
    ))
  # wait a tinybit before returning link (backoff)
  time.sleep(1)
  result = resp.json()
  return dict(result, link=enrichr_link + '/enrich?dataset=' + resp.json()['shortId'])

def enrichr_get_top_results(userListId, bg, enrichr_link='https://amp.pharm.mssm.edu/Enrichr'):
  import time, requests
  time.sleep(1)
  resp = requests.get(enrichr_link + '/enrich?userListId={}&backgroundType={}'.format(userListId, bg))
  if resp.status_code != 200:
    raise Exception('Enrichr failed with status {}: {}'.format(
      resp.status_code,
      resp.text,
    ))
  time.sleep(1)
  return pd.DataFrame(resp.json()[bg], columns=['rank', 'term', 'pvalue', 'zscore', 'combinedscore', 'overlapping_genes', 'adjusted_pvalue', '', ''])


# Merge data
df_clustered_tsne = pd.merge(left=df_clusters, left_on='Barcode', right=df_tsne, right_on='Barcode')
df_clustered_pca = pd.merge(left=df_clusters, left_on='Barcode', right=df_pca, right_on='Barcode')

# Grab ncbi symbols
ncbi = pd.read_csv('ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz', sep='\t')
# Ensure nulls are treated as such
ncbi = ncbi.applymap(lambda v: float('nan') if type(v) == str and v == '-' else v)
# Break up lists
split_list = lambda v: v.split('|') if type(v) == str else []
ncbi['dbXrefs'] = ncbi['dbXrefs'].apply(split_list)
ncbi['Synonyms'] = ncbi['Synonyms'].apply(split_list)
ncbi['LocusTag'] = ncbi['LocusTag'].apply(split_list)
ncbi['Other_designations'] = ncbi['Other_designations'].apply(split_list)
# Map existing entities to NCBI Genes
ncbi_lookup = {
  sym.upper(): row['Symbol'].upper()
  for _, row in ncbi.iterrows()
  for sym in [row['Symbol']] + row['Synonyms']
}
df['Symbol'] = df['Feature Name'].map(lambda s: ncbi_lookup.get(s.upper()))

# Get top Genes for each cluster
top_genes = {}
for cluster in map(str, range(1, 10+1)):
  fc_col = 'Cluster %s Log2 fold change' % (cluster)
  p_col = 'Cluster %s Adjusted p value' % (cluster)
  # significant and positive fold change sorted by p value
  up_genes = df.loc[
    df[((df[p_col] < 0.05) & (df[fc_col] > 0))][p_col].sort_values().index,
    'Symbol'
  ].iloc[:n_genes].dropna().values
  # significant and negative fold change sorted by p value
  dn_genes = df.loc[
    df[((df[p_col] < 0.05) & (df[fc_col] < 0))][p_col].sort_values().index,
    'Symbol'
  ].iloc[:n_genes].dropna().values
  # save results
  top_genes[cluster] = (up_genes, dn_genes)


# Get Enrichr links for each cluster
enrichr_links = {}

for cluster, (up_genes, dn_genes) in top_genes.items():
  up_link, dn_link = None, None
  if up_genes.size:
    up_link = enrichr_link_from_genes(up_genes, 'cluster %s up' % (cluster))
    # display_link_inline(up_link['link'])
  else:
    print('cluster %s up: empty' % (cluster))
  if dn_genes.size:
    dn_link = enrichr_link_from_genes(dn_genes, 'cluster %s down' % (cluster))
    # display_link_inline(dn_link['link'])
  else:
    print('cluster %s down: empty' % (cluster))
  enrichr_links[cluster] = (up_link, dn_link)

# Grab top results for each cluster
all_results = []
for cluster, (up_link, dn_link) in enrichr_links.items():
  for link_type, link in [('up', up_link), ('down', dn_link)]:
    if link is None:
      continue
    for category, libraries in useful_libs.items():
      for library in libraries:
        try:
          results = enrichr_get_top_results(link['userListId'], library).sort_values('pvalue').iloc[:top_n_results]
          results['link'] = link['link']
          results['library'] = library
          results['category'] = category
          results['direction'] = link_type
          results['cluster'] = cluster
          all_results.append(results)
        except:
          print('{}: {} {} {} cluster {} failed, continuing'.format(link, library, category, direction, cluster))

df_all_results = pd.concat(all_results)

df.to_csv('data/df.tsv', sep='\t')
df_clustered_tsne.to_csv('data/df_tsne.tsv', sep='\t')
df_all_results.to_csv('data/df_enrich.tsv', sep='\t')
