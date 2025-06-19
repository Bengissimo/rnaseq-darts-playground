## Time series forecasting of gene expressions

Forecast gene expression trends over time using Darts and RNA-seq data.

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

### Dataset: GSE253864 - Human Blood Transcriptome During Bed Rest

This repository includes reprocessed data derived from the publicly available dataset [GSE253864](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253864), hosted on [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/).

#### Citation

Archer SN, Möller-Levet C, Bonmatí-Carrión MÁ, Laing EE et al.  
**Extensive dynamic changes in the human transcriptome and its circadian organization during prolonged bed rest.**  
_iScience._ 2024 Mar 15;27(3):109331.  
PMID: [38487016](https://pubmed.ncbi.nlm.nih.gov/38487016/)  
DOI: [10.1016/j.isci.2024.109331](https://doi.org/10.1016/j.isci.2024.109331)
