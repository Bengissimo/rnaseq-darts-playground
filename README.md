## Time series forecasting of gene expressions

I built this project to practice time series forecasting using the Darts Python library, applying it to real biological data from a public RNA-seq microarray dataset (NCBI GEO). I wanted to explore how forecasting models perform on gene expression profiles that change over time.

To test the models, I focused primarily on PER1, a well-known circadian clock gene with a 24-hour cycling expression pattern. Its predictable rhythm made it a useful benchmark for evaluating different forecasting methods.

While the project began as a learning exercise, it can be applied to any gene of interest to visualize expression trends andexperiment with forecasting techniques on real omics data.

## How to run

### 1. **Install dependencies**

Make sure you have Python 3.8+ and install required packages:

```bash
pip install -r requirements.txt
```

### 2. **Run the forecasting**

Use the command line to run forecasting on selected genes:

```bash
python3 time_series.py --genes GENE1 GENE2 --mode classical --plots
```

* `--genes`: space-separated list of gene symbols (e.g. `PER1 NR1D1`)
* `--mode`: forecasting mode: `classical` (default) or `deep_learning`
* `--plots`: (optional) show forecast plots

### Example:

```bash
python3 time_series.py --genes PER1 --mode deep_learning --plots
```

## Dataset

**GSE253864 - Human Blood Transcriptome During Bed Rest**
This repository includes reprocessed data derived from the publicly available dataset [GSE253864](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253864), hosted on [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/).

### Data Preparation
To generate a continuous time series for forecasting, I prepared gene expression data from microarray measurements by grouping and aggregating values across samples. Specifically, I computed the median expression values per time point for each gene across different biological groups and experimental replicates. Using known sample metadata and starting dates, I reconstructed a unified timeline covering approximately 60 days, assigning each expression point a corresponding timestamp. This allowed me to convert the data into a clean, datetime-indexed format suitable for time series modeling with Darts.

### Citation

Archer SN, Möller-Levet C, Bonmatí-Carrión MÁ, Laing EE et al.  
**Extensive dynamic changes in the human transcriptome and its circadian organization during prolonged bed rest.**  
_iScience._ 2024 Mar 15;27(3):109331.  
PMID: [38487016](https://pubmed.ncbi.nlm.nih.gov/38487016/)  
DOI: [10.1016/j.isci.2024.109331](https://doi.org/10.1016/j.isci.2024.109331)
