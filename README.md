# Single-cell RNA-seq analysis code — MOG antibody-associated encephalitis

Analysis code accompanying the manuscript:

> **<Manuscript title — please fill in>**
> <Author list>. Manuscript ID: RX2026-04076-W.

This repository contains the R / Python / Perl scripts used for quality control,
integration, clustering, cell-type annotation, differential expression, cell–cell
communication, transcription-factor activity and pseudotime analyses of
single-cell RNA-seq data from patients with MOG antibody-associated
encephalitis (MOGAD) and healthy controls (HC).

Raw sequencing data and per-cell expression matrices are **not** included here
(see *Data availability*).

---

## Repository structure

| Directory | Script | Purpose |
|---|---|---|
| `基础/` | `sc_preprocessing.r` | Quality control, doublet removal (DoubletFinder), ambient RNA correction (SoupX), integration (Harmony), clustering and UMAP |
| `基础/` | `sc_cells_analysis.r` | Cell-type annotation, marker detection, differential abundance and DEG analysis across cell types |
| `基础/` | `sc_genes_analysis.r` | Gene-level analyses, including pseudobulk differential expression between MOG and HC |
| `基础/` | `sc_cell_pct_plot.r` | Cell-composition plotting |
| `基础/` | `sc_correlation.r` | Correlation analyses across samples / cell types |
| `基础/` | `scRNAtoolVis_plot.r` | Volcano and summary plots (scRNAtoolVis) |
| `基础/` | `kegg_go_pipe_one.pl` | GO / KEGG functional enrichment (Perl pipeline wrapper) |
| `拟时序/` | `monocle2_basic_v2.r` | Pseudotime / trajectory analysis with Monocle 2 |
| `通讯/` | `sc_cellchat_bycelltype_1_sample.r` | Cell–Cell communication (CellChat) within a single sample, per cell type |
| `通讯/` | `sc_cellchat_bycelltype_2_compare.r` | Comparative CellChat analysis between MOG and HC |
| `转录因子/` | `Seurat2loom.py` | Convert Seurat objects to `.loom` for SCENIC |
| `转录因子/` | `SCENIC_Plot_bycelltype_v2.r` | SCENIC regulon activity, per cell type |
| `转录因子/` | `AUC_plot_v2.r` | AUCell regulon-activity plots |
| `转录因子/` | `sc_genes_plot.r` | Expression plots for selected genes / regulons |

---

## Requirements

- **R** ≥ 4.1 (scripts were developed and tested on R 4.x)
- **Python** ≥ 3.8 (only for `转录因子/Seurat2loom.py`)
- **Perl** ≥ 5 (only for `基础/kegg_go_pipe_one.pl`)

### R packages

Core single-cell:

```
Seurat, SeuratDisk, SeuratWrappers, SingleCellExperiment, scater, scRNAseq,
harmony, DoubletFinder, SoupX, clustree, presto, scRNAtoolVis
```

Downstream analyses:

```
monocle, CellChat, SCENIC, SCopeLoomR, AUCell, NMF, hdf5r
```

Statistics and visualisation:

```
limma, metap, BiocParallel, doParallel, future.apply, data.table, tidyverse,
dplyr, stringr, reshape2, scales, gtools, grid, grDevices, stats, KernSmooth,
ggplot2, ggpubr, ggrepel, ggsci, ggthemes, ggalluvial, gplots, cowplot,
patchwork, pheatmap, ComplexHeatmap, circlize, RColorBrewer, plotly, optparse
```

### Python packages

```
scanpy, loompy, numpy, optparse
```

---

## How to run

Scripts are numbered by analysis stage and are designed to be run in the
following order. Paths to input objects are passed as command-line arguments
(`optparse`); edit the `setwd()` / input-path variables at the top of each
script to match your local environment.

```bash
# 1. Preprocessing, integration and clustering
Rscript 基础/sc_preprocessing.r

# 2. Cell-type annotation and differential analyses
Rscript 基础/sc_cells_analysis.r
Rscript 基础/sc_genes_analysis.r

# 3. Descriptive / summary figures
Rscript 基础/sc_cell_pct_plot.r
Rscript 基础/sc_correlation.r
Rscript 基础/scRNAtoolVis_plot.r

# 4. Pseudotime
Rscript 拟时序/monocle2_basic_v2.r

# 5. Cell-cell communication
Rscript 通讯/sc_cellchat_bycelltype_1_sample.r
Rscript 通讯/sc_cellchat_bycelltype_2_compare.r

# 6. Transcription-factor (SCENIC) analysis
python 转录因子/Seurat2loom.py
Rscript 转录因子/SCENIC_Plot_bycelltype_v2.r
Rscript 转录因子/AUC_plot_v2.r
Rscript 转录因子/sc_genes_plot.r
```

Intermediate `.rds` / `.loom` objects are written to the working directory and
are not tracked in this repository.

---

## Data availability

The single-cell RNA-seq data supporting this study have been deposited in
**<repository, e.g. GEO / GSA / SRA>** under accession **<accession ID>**
[add URL once available].

Because the data derive from human participants, individual-level expression
matrices are not distributed through this repository. Processed data may be
obtained from the corresponding author upon reasonable request, subject to
institutional ethics approval.

---

## License

<Choose one — e.g. MIT / GPL-3.0 / CC BY 4.0>. See `LICENSE`.

---

## Contact

Corresponding author: **<name>** — <email>
Code maintained by: **Tong Qianwen** — qianwentong@rxbio.cc
