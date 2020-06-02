#%%
import os
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
from IPython.display import display
from matplotlib import pyplot as plt
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from maayanlab_bioinformatics.dge import characteristic_direction
from maayanlab_bioinformatics.normalization import log2_normalize, filter_by_var, zscore_normalize
from maayanlab_bioinformatics.utils import merge

#%%
df_data = pd.read_csv('data.csv', index_col=0)
df_data.columns = df_data.columns.astype(str)
df_data.columns.name = 'Barcode'
metadata = pd.read_csv('metadata.csv', index_col=0)
metadata.index = metadata.index.astype(str)
metadata['Labels'] = metadata['Labels'].astype(str)

#%%
df_library_size = pd.DataFrame(
    {
        'n_reads': df_data[df_data > 0].count(),
        'log_n_reads': np.log2(df_data[df_data > 0].count() + 1),
        'n_expressed_genes': df_data.sum(),
    }
).sort_values('n_reads', ascending=False)

display(df_library_size.head())
sns.distplot(df_data.iloc[0, :]); plt.show()
sns.distplot(df_data.iloc[:, 0]); plt.show()

#%%
df_data_norm = filter_by_var(df_data)
df_data_norm = log2_normalize(df_data_norm)
df_data_norm = zscore_normalize(df_data_norm)

#%%
sns.distplot(df_data_norm.iloc[0, :]); plt.show()
sns.distplot(df_data_norm.iloc[:, 0]); plt.show()

#%%
data_norm_pca = PCA(
  random_state=42,
)
data_norm_pca.fit(df_data_norm.values.T)
df_data_norm_pca = pd.DataFrame(
    data_norm_pca.transform(df_data_norm.values.T),
    index=df_data_norm.T.index
)
df_data_norm_pca.columns = [
    f'PCA-{c}' # ({r:.3f})'
    for c, r in zip(df_data_norm_pca.columns, data_norm_pca.explained_variance_ratio_)
]

#%%
px.scatter(
  merge(
    df_data_norm_pca,
    df_library_size,
    metadata,
  ),
  x=df_data_norm_pca.columns[0],
  y=df_data_norm_pca.columns[1],
  size='n_reads',
  size_max=8,
  color='Labels',
  hover_data=[df_data_norm.columns],
)

#%%
data_norm_umap = UMAP(
  random_state=42,
  n_components=2,
  n_neighbors=30,
  metric='cosine',
  min_dist=0.3,
)
data_norm_umap.fit(df_data_norm_pca.iloc[:, :10].values)
df_data_norm_umap = pd.DataFrame(
  data_norm_umap.transform(df_data_norm_pca.iloc[:, :10].values),
  columns=['UMAP-1', 'UMAP-2'],
  index=df_data_norm_pca.index,
)

#%%
px.scatter(
  merge(
    df_data_norm_umap,
    df_library_size,
    metadata,
  ),
  x=df_data_norm_umap.columns[0],
  y=df_data_norm_umap.columns[1],
  size='n_reads',
  size_max=8,
  color='Labels',
  hover_data=[df_data_norm.columns],
)

#%%
silhouette_scores = {}
for n in range(2, 25):
    y_pred = KMeans(n_clusters=n, random_state=42).fit_predict(df_data_norm_umap.values)
    silhouette_scores[n] = silhouette_score(df_data_norm_umap.values, y_pred, metric='cosine')

silhouette_scores = pd.DataFrame([
    {'N Clusters': k, 'Silhouette Score': v}
    for k, v in silhouette_scores.items()
])
best = silhouette_scores.sort_values('Silhouette Score').iloc[-1]
silhouette_scores

#%%
plt.plot(silhouette_scores['N Clusters'], silhouette_scores['Silhouette Score'])
plt.scatter([best['N Clusters']], [best['Silhouette Score']], label='Best')
plt.legend()
plt.title('Cluster size selection')
plt.ylabel('Silhouette Score')
plt.xlabel('Number of Clusters')
plt.show()

#%%
km = KMeans(n_clusters=int(best['N Clusters']), random_state=42)
df_data_norm_km = pd.DataFrame({
    'Cluster': [
        str(c)
        for c in km.fit_predict(df_data_norm_umap.values)
    ]
}, index=df_data_norm_umap.index)

#%%
# Perform differential expression for each cluter
diff_expr = {}
for cluster, samples in df_data_norm_km.groupby('Cluster'):
  diff_expr[cluster] = characteristic_direction(
    # expression outside of this cluster
    df_data_norm.loc[:, df_data_norm.columns.difference(samples.index)],
    # expression in this cluster
    df_data_norm.loc[:, samples.index],
  )['CD-coefficient']

df_diff_expr = pd.DataFrame(diff_expr)
df_diff_expr.index.name = 'Feature Name'

#%%
df_diff_expr['Cluster 0'].sort_values(ascending=True)
#%%
aucs = {}
for cluster, samples in df_data_norm_km.groupby('Cluster'):
  aucs[cluster] = {}
  for feature in metadata.columns:
    lr = LogisticRegression()
    X = metadata.loc[:, feature].values[:, np.newaxis]
    y_true = (df_data_norm_km['Cluster'] == cluster).astype(np.int64)
    lr.fit(X, y_true)
    y_score = lr.predict_proba(X)[:, 1]
    aucs[f"Cluster {cluster} CD"][feature] = roc_auc_score(y_true, y_score)

pd_aucs = pd.DataFrame(aucs)

#%%
# /clustering/graphclust/clusters.csv
#   Barcode,Cluster
os.makedirs('clustering/graphclust', exist_ok=True)
df_data_norm_km.to_csv('clustering/graphclust/clusters.csv')
# /diffexp/graphclust/differential_expression.csv
#   Feature Name,Cluster 2 Log2 fold change,Cluster 0 Log2 fold change,Cluster 1 Log2 fold change,Cluster 3 Log2 fold change,Cluster 2 Adjusted p value,Cluster 0 Adjusted p value,Cluster 1 Adjusted p value,Cluster 3 Adjusted p value,Cluster 2 Mean Counts,Cluster 0 Mean Counts,Cluster 1 Mean Counts,Cluster 3 Mean Counts,Feature Name
os.makedirs('diffexp/graphclust', exist_ok=True)
df_diff_expr.to_csv('diffexp/graphclust/differential_expression.csv')
# /pca/10_components/projection.csv
#   Barcode,PC-1,PC-2,PC-3,PC-4,PC-5,PC-6,PC-7,PC-8,PC-9,PC-10
os.makedirs('pca/10_components', exist_ok=True)
df_data_norm_pca.to_csv('pca/10_components/projection.csv')
# /umap/2_components/projection.csv
#   Barcode,UMAP-1,UMAP-2
os.makedirs('umap/2_components', exist_ok=True)
df_data_norm_umap.to_csv('umap/2_components/projection.csv')
# /cluster_aucs.csv
#   ,2,0,1,3
pd_aucs.to_csv('cluster_aucs.csv')
