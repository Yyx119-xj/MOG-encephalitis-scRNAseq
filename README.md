# Single-cell transcriptome reveals immune cell subset dysregulation and intercellular signaling remodeling in MOG encephalitis

Analysis code accompanying the manuscript:

> **Single-cell transcriptome reveals immune cell subset dysregulation and
> intercellular signaling remodeling in MOG encephalitis**
>
> Yuxuan Yong<sup>1,#</sup>, Ting Liu<sup>2,#</sup>, Mengdie Hao<sup>1</sup>,
> Feifei Pei<sup>1</sup>, Xian Su<sup>3</sup>, Ye Du<sup>1</sup>,
> Yixuan Ma<sup>1</sup>, Fangfang Bi<sup>4,\*</sup>, Xinling Yang<sup>1,\*</sup>
>
> <sup>1</sup>Department of Neurology, Second Affiliated Hospital of Xinjiang Medical University, Urumqi, Xinjiang, China
>
> <sup>2</sup>Department of Neurology, Affiliated Traditional Chinese Medicine Hospital of Xinjiang Medical University, Urumqi, Xinjiang, China
>
> <sup>3</sup>Department of Neurology, The 474 Hospital of Chinese PLA, Urumqi, Xinjiang, China
>
> <sup>4</sup>Department of Neurology, The Fifth Affiliated Hospital of Sun Yat-sen University, Zhuhai, Guangdong, China
>
> <sup>#</sup>These authors contributed equally to this work.
>
> <sup>\*</sup>Co-corresponding authors: Fangfang Bi (biff3@mail.sysu.edu.cn);
> Xinling Yang (yangxinling2014@163.com)

This repository contains the R / Python / Perl scripts used for quality control,
integration, clustering, cell-type annotation, differential expression, cell–cell
communication, transcription-factor activity and pseudotime analyses of
single-cell RNA-seq data from patients with MOG antibody-associated
encephalitis (MOG encephalitis) and healthy controls (HC).

Raw sequencing data and per-cell expression matrices are **not** included here
(see [Data availability](#data-availability)).

---

## Repository structure

| Directory | Script | Purpose |
|---|---|---|
| `basic/` | `sc_preprocessing.r` | Quality control, doublet removal (DoubletFinder), ambient RNA correction (SoupX), integration (Harmony), clustering and UMAP |
| `basic/` | `sc_cells_analysis.r` | Cell-type annotation, marker detection, differential abundance and DEG analysis across cell types |
| `basic/` | `sc_genes_analysis.r` | Gene-level analyses, including pseudobulk differential expression between MOG and HC |
| `basic/` | `sc_cell_pct_plot.r` | Cell-composition plotting |
| `basic/` | `sc_correlation.r` | Correlation analyses across samples / cell types |
| `basic/` | `scRNAtoolVis_plot.r` | Volcano and summary plots (scRNAtoolVis) |
| `basic/` | `kegg_go_pipe_one.pl` | GO / KEGG functional enrichment (Perl pipeline wrapper) |
| `monocle2/` | `monocle2_basic_v2.r` | Pseudotime / trajectory analysis with Monocle 2 |
| `cellchat/` | `sc_cellchat_bycelltype_1_sample.r` | Cell–cell communication (CellChat) within a single sample, per cell type |
| `cellchat/` | `sc_cellchat_bycelltype_2_compare.r` | Comparative CellChat analysis between MOG and HC |
| `SCENIC/` | `Seurat2loom.py` | Convert Seurat objects to `.loom` for SCENIC |
| `SCENIC/` | `SCENIC_Plot_bycelltype_v2.r` | SCENIC regulon activity, per cell type |
| `SCENIC/` | `AUC_plot_v2.r` | AUCell regulon-activity plots |
| `SCENIC/` | `sc_genes_plot.r` | Expression plots for selected genes / regulons |

---

## Requirements

- **R** ≥ 4.1 (scripts were developed and tested on R 4.x)
- **Python** ≥ 3.8 (only for `SCENIC/Seurat2loom.py`)
- **Perl** ≥ 5 (only for `basic/kegg_go_pipe_one.pl`)

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

Scripts are organised by analysis stage and are designed to be run in the
following order. Paths to input objects are passed as command-line arguments
(`optparse`); edit the `setwd()` / input-path variables at the top of each
script to match your local environment.

```bash
# 1. Preprocessing, integration and clustering
Rscript basic/sc_preprocessing.r

# 2. Cell-type annotation and differential analyses
Rscript basic/sc_cells_analysis.r
Rscript basic/sc_genes_analysis.r

# 3. Descriptive / summary figures
Rscript basic/sc_cell_pct_plot.r
Rscript basic/sc_correlation.r
Rscript basic/scRNAtoolVis_plot.r

# 4. Pseudotime
Rscript monocle2/monocle2_basic_v2.r

# 5. Cell-cell communication
Rscript cellchat/sc_cellchat_bycelltype_1_sample.r
Rscript cellchat/sc_cellchat_bycelltype_2_compare.r

# 6. Transcription-factor (SCENIC) analysis
python SCENIC/Seurat2loom.py
Rscript SCENIC/SCENIC_Plot_bycelltype_v2.r
Rscript SCENIC/AUC_plot_v2.r
Rscript SCENIC/sc_genes_plot.r
```

Intermediate `.rds` / `.loom` objects are written to the working directory and
are not tracked in this repository.

---

## Data availability

The single-cell RNA-seq data supporting this study have been deposited in the
**NCBI Gene Expression Omnibus (GEO)** under accession
**[GSE353303](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353303)**.

The series is private during peer review and will be released upon publication.
<!-- TODO: insert the GEO reviewer token link here before submission, e.g.
     https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353303&token=xxxxxxxx -->

Because the data derive from human participants, individual-level expression
matrices are not distributed through this repository. Processed data may be
obtained from the corresponding author upon reasonable request, subject to
institutional ethics approval.

---

## Citation

If you use this code, please cite the manuscript above. A permanent DOI for
this repository will be added here once the paper is accepted.

---

## License

This code is released under the **MIT License** — see [LICENSE](LICENSE).
You are free to use, modify and redistribute it, provided that the copyright
notice is retained.

---

## Contact

Co-corresponding authors:

- **Fangfang Bi** — biff3@mail.sysu.edu.cn
- **Xinling Yang** — yangxinling2014@163.com

For questions about the analysis, please open an issue in this repository.
