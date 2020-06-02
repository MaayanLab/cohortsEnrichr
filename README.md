# cohortsEnrichr

A Plotly Dash application for interactive exploration of 10x Genomics Output

This enables 10X Genomics cluster exploration by means of interactive cluster enrichment analysis. The top enrichment results are pre-computed with `init.py` prior to using the interactive web application `app.py`.

## Preparation
### Setup
```bash
python3 -m venv venv
source venv/bin/activate
pip3 install -r requirements.txt
```

### Prepare data files
```bash
source venv/bin/activate
mkdir -p data
python3 init.py 10x_output_directory/outs/analysis data
```

## Usage
### Run application locally
```bash
source venv/bin/activate
DATA=data python3 app.py
```

### Run application with docker-compose
```bash
DATA=data docker-compose up
```

## Data
The input data format for init should be of the form:
```
/clustering/graphclust/clusters.csv
  Barcode,Cluster
/diffexp/graphclust/differential_expression.csv
  Feature ID,Cluster 2 Log2 fold change,Cluster 0 Log2 fold change,Cluster 1 Log2 fold change,Cluster 3 Log2 fold change,Cluster 2 Adjusted p value,Cluster 0 Adjusted p lue,Cluster 1 Adjusted p value,Cluster 3 Adjusted p value,Cluster 2 Mean Counts,Cluster 0 Mean Counts,Cluster 1 Mean Counts,Cluster 3 Mean Counts,Feature Name
/pca/10_components/projection.csv
  Barcode,PC-1,PC-2,PC-3,PC-4,PC-5,PC-6,PC-7,PC-8,PC-9,PC-10
/tsne/2_components/projection.csv
  Barcode,TSNE-1,TSNE-2
/cluster_aucs.csv
  ,2,0,1,3
/metadata.csv
  Barcode,type_subject,...
```

An example of turning a `data.csv` and `metadata.csv` into this output is available in `example/` as well as the resulting files such that `example/` can also be used with `init.py` for `app.py` with `DATA=example python3 app.py`.
