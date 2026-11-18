# RNA-seq QC and Exploratory Analysis in Python

A reusable command-line pipeline for quality checks and exploratory analysis of **raw RNA-seq gene count matrices**. Designed for researchers who want reproducible tables and plots without editing the analysis functions.

> **Scope:** This tool starts from a gene-by-sample count matrix. It does not process FASTQ files, perform read alignment, estimate differential expression, or correct batch effects.

## Features

- Validate nonnegative, integer count data; report missing values and duplicate identifiers
- Summarize library sizes and detected genes per sample
- Filter genes by minimum count in a minimum number of samples
- Calculate counts per million (CPM) using **original, unfiltered library totals**, then log2(CPM + 1)
- Plot sample-wise expression distributions
- Compute Pearson or Spearman sample correlations and heatmaps
- Perform exploratory PCA, optionally colored by metadata-defined groups
- Export tables as CSV and figures as PNG

## Installation

Python 3.10+ is recommended. In a terminal inside the project directory:

```bash
python -m pip install -r requirements.txt
```

## Quick start

Run the synthetic example:

```bash
python pipeline_rnaseq.py --rna examples/example_counts.csv --metadata examples/example_metadata.csv --coluna-amostra sample_id --coluna-grupo group --saida example_results --minimo-amostras 2
```

Without metadata:

```bash
python pipeline_rnaseq.py --rna examples/example_counts.csv --saida example_results_no_metadata --minimo-amostras 2
```

To view all options:

```bash
python pipeline_rnaseq.py --help
```

## Input files

**Count matrix:** genes in rows, samples in columns, gene identifiers in the first column, and raw nonnegative integer counts in all other columns. Supported: `.csv` (comma-separated), `.tsv` / `.tab` (tab-separated), and whitespace-delimited `.txt`. Whitespace-delimited TXT files with a header missing the gene-ID column name are automatically detected when the data row has one extra field.

**Optional metadata:** CSV with one sample-ID column and one group column. Column names are configurable with `--coluna-amostra` and `--coluna-grupo`. Sample IDs must match count-matrix column names. Additional metadata samples are allowed; missing group labels for count-matrix samples are marked as ungrouped.

## Command-line settings

| Option | Default | Meaning |
|---|---|---|
| `--rna` | `rawcountAbate.txt` | Input count matrix |
| `--metadata` | none | Optional sample metadata CSV |
| `--saida` | `Resultados_RNAseq` | Output directory |
| `--separador` | `auto` | `auto`, `virgula`, `tab`, or `espacos` |
| `--cabecalho-sem-id` | off | Explicitly indicate unnamed gene-ID header for whitespace TXT |
| `--contagem-minima` | `10` | Minimum raw count |
| `--minimo-amostras` | `3` | Number of samples meeting minimum count |
| `--coluna-amostra` | `#NAME` | Metadata sample ID column |
| `--coluna-grupo` | `class` | Metadata group column |
| `--correlacao` | `pearson` | `pearson` or `spearman` |

The default settings can also be edited near the top of `pipeline_rnaseq.py`. Command-line options override those defaults. For datasets with fewer than three samples, PCA is skipped; other outputs may still be produced.

## Output files

The output folder contains `QC_validacao.csv`, `QC_resumo_amostras.csv`, `QC_tamanho_bibliotecas.png`, `resumo_processamento.csv`, `RNAseq_log2_CPM.csv`, `boxplot_log2_CPM.png`, `matriz_correlacao.csv`, `correlacao_media.csv`, `heatmap_correlacao_0a1.png`, `heatmap_correlacao_ajustado.png`, `PCA_coordenadas.csv`, `PCA_variancia_explicada.csv`, and `PCA_RNAseq.png` (or `PCA_RNAseq_por_grupo.png` when metadata are supplied). When metadata are supplied, `PCA_coordenadas_com_grupos.csv` is also written.

## Notes on interpretation

CPM and log2(CPM + 1) are useful for basic exploratory plots, but they are **not a replacement for differential-expression normalization and modeling**. PCA is computed from the centered log2-CPM matrix without scaling each gene to unit variance. The filtering rule and plots are exploratory choices and may require adaptation for particular study designs.

## Data sharing and attribution

The `examples/` files are synthetic and contain no real study data. Do not publish private, restricted, or unpublished experimental datasets without permission. This repository intentionally contains no institutional affiliation or author attribution; add those details after confirming the appropriate credits. No open-source license is included yet; choose one with the relevant collaborators or institution before inviting reuse or redistribution.

## Development

This is an initial reusable release. To contribute or report issues, use the GitHub repository's Issues section after publication.
