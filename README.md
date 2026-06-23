# scRNA-seq Analysis of Immune Cell Transcriptomics Under Microgravity

Single-cell RNA sequencing (scRNA-seq) analysis of human peripheral blood mononuclear cells (PBMCs) exposed to simulated microgravity (µG) compared to normal gravity (1G) controls, with and without immune stimulation.

**Author:** Anastasiya Pautarak  
**Institution:** Jagiellonian University (UJ), Faculty of Biochemistry, Biophysics and Biotechnology, 2026  
**Supervisor:** dr Adrian Kania  
**GEO Accession:** [GSE218936](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218936)

---

## Project Overview

This study investigates how microgravity conditions affect the transcriptomic landscape of human immune cells. PBMCs from two donors were cultured under four experimental conditions:

| Condition | Gravity | Stimulation |
|-----------|---------|-------------|
| `1G_stim` | Normal (1G) | PHA-stimulated |
| `1G_unstim` | Normal (1G) | Unstimulated |
| `uG_stim` | Microgravity (µG) | PHA-stimulated |
| `uG_unstim` | Microgravity (µG) | Unstimulated |

**Key analyses performed:**
- Quality control, filtering and normalization of raw scRNA-seq counts
- Automated cell type annotation using CellTypist (`Immune_All_Low` model)
- UMAP dimensionality reduction and graph-based clustering (Leiden / Louvain)
- Pseudobulk PCA for sample-level correlation
- Differential gene expression (Wilcoxon test, Volcano plots)
- Gene Ontology enrichment analysis (GSEApy / Enrichr)
- Jaccard index-based transcriptomic variability per cell type

---

## Repository Structure

```
.
├── 1_data_and_filtration.ipynb   # QC, normalization, PCA, cell type annotation
├── 2_umap_analisys.ipynb         # UMAP, clustering, batch effect, pseudobulk PCA
├── 3_expression_analysis.ipynb   # DGE, GO, Jaccard, stress score, PAGA, LIANA
├── env_scRNA.yml                 # Conda environment specification
├── figures/                      # Output figures (PNG)
└── README.md
```

> **Note:** Raw data (`.h5ad` checkpoints, raw count matrices) are not included in this repository due to file size. See [Data](#data) section below.

---

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- ~32 GB RAM recommended (dataset: ~60 000 cells × 36 601 genes)
- Python 3.11

---

## Environment Setup

### 1. Clone the repository

```bash
git clone https://github.com/anastasiyapautarak/scRNA_microG_immunity.git
cd scRNA_microG_immunity
```

### 2. Create and activate the conda environment

```bash
conda env create -f env_scRNA.yml
conda activate scRNA_microG_immunity
```

This will install all required packages including:
- `scanpy`, `anndata` — core scRNA-seq analysis
- `celltypist` — automated cell type annotation
- `liana` — cell-cell communication
- `gseapy` — Gene Ontology enrichment
- `umap-learn`, `leidenalg`, `louvain` — dimensionality reduction and clustering
- `jupyterlab` — notebook interface

### 3. Launch JupyterLab

```bash
jupyter lab
```

---

## Data

Raw data is available from NCBI GEO under accession **[GSE218936](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218936)**.

Download the 8 sample folders (10x Genomics Cell Ranger output format) and place them in the project directory with the following naming convention:

```
SC_040522_1G_stim/
SC_040522_1G_unstim/
SC_040522_uG_stim/
SC_040522_uG_unstim/
SC_050322_1G_stim/
SC_050322_1G_unstim/
SC_050322_uG_stim/
SC_050322_uG_unstim/
```

Each folder must contain the standard Cell Ranger output files: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`.

---

## Reproducing the Analysis

Run the notebooks **in order**. Each notebook saves a checkpoint (`.h5ad` file) that is loaded by the next one.

### Notebook 1 — `1_data_and_filtration.ipynb`

Loads raw count matrices from all 8 samples, parses experimental metadata from folder names, and performs:

- **QC filtering:** min 200 genes/cell, max 4000 genes/cell, <10% mitochondrial reads, <5% hemoglobin genes
- **Normalization:** total count normalization to 10 000, log1p transformation
- **Feature selection:** highly variable genes (HVGs)
- **PCA** and elbow plot for component selection
- **Cell type annotation** via CellTypist (`Immune_All_Low` model)

Saves checkpoint: `data_checkpoints/adata_po_QC_i_PCA.h5ad`

### Notebook 2 — `2_umap_analisys.ipynb`

Loads the QC checkpoint and computes:

- **KNN graph** and **UMAP** embedding
- **Clustering** with Leiden and Louvain algorithms at multiple resolutions
- **Split UMAP** per sample to assess batch effects
- **Pseudobulk PCA** with Pearson correlation heatmap between samples

Saves checkpoint: `data_checkpoints/dane_po_anotacji_z_umap.h5ad`

### Notebook 3 — `3_expression_analysis.ipynb`

Loads the annotated UMAP checkpoint and runs:

- **Mean/median expression statistics** per condition and cell type → `figures/wyniki_ekspresji/`
- **Gene Ontology enrichment** (top 100 differentially expressed genes, Enrichr API) → `figures/gene_ontology/`
- **Jaccard index** of active gene sets between µG and 1G per cell type
- **Stress gene scoring** (HIF1A, HSP90AA1, HSPA1A, SOD2, etc.) with UMAP and violin plots
- **Volcano plot** (Wilcoxon test, µG vs 1G, FDR < 0.05, |log2FC| > 1)
- **DGE balance** (up/downregulated genes per cell type in µG)

> ⚠️ Notebook 3 is computationally intensive. Recommended minimum: 32 GB RAM, 8 CPU cores.

---

## Output Figures

| File | Description |
|------|-------------|
| `figures/elbow_plot_PCA.png` | PCA elbow plot for component selection |
| `figures/8_probek_UMAP.png` | UMAP split by sample |
| `figures/UMAP_glowny_cell_type.png` | Main UMAP coloured by cell type |
| `figures/UMAP_bez_labeli.png` | UMAP without labels |
| `figures/porownanie_rozdzielczosci_leiden.png` | Leiden clustering at multiple resolutions |
| `figures/porownanie_rozdzielczosci_louvain.png` | Louvain clustering at multiple resolutions |
| `figures/korelacja_stim.png` | Pseudobulk PCA correlation — stimulated |
| `figures/korelacja_unstim.png` | Pseudobulk PCA correlation — unstimulated |
| `figures/1_UMAP_stress_score.png` | UMAP coloured by cellular stress score |
| `figures/1b_violin_stress_score.png` | Stress score per experimental condition |
| `figures/2_PAGA_network.png` | PAGA connectivity graph |
| `figures/3_UMAP_cell_types_paga.png` | UMAP with PAGA-initialized layout |
| `figures/3_UMAP_pseudotime.png` | Diffusion pseudotime on UMAP |
| `figures/4_Volcano_uG_vs_1G.png` | Volcano plot: µG vs 1G |
| `figures/zmiennosc_komorek_jaccard.csv` | Jaccard variability scores per cell type |
