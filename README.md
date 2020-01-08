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
python3 init.py 10x_output_directory/outs/analysis
```

## Usage
### Run application locally
```bash
source venv/bin/activate
python3 app.py
```

### Run application with docker-compose
```bash
docker-compose up
```
