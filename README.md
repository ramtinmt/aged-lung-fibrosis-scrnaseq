# Nerandomilast in the aged fibrotic lung — single-cell RNA-seq

Single-cell analysis of the *Sftpc*<sup>I73T</sup> mouse model of pulmonary
fibrosis, testing whether the PDE4B inhibitor **nerandomilast** resolves
fibrosis differently in **young (3-month)** versus **aged (18-month)** lungs.

<!-- TODO: one sentence on what you actually found. Write this last. -->

## The question

Nerandomilast is in trials for IPF, a disease of the aging lung — but efficacy
is usually established in young animals. This experiment asks whether the drug's
effect on lung cell composition and cell–cell signalling is **age-dependent**.

## Experimental design

8 samples: 3 conditions × 2 ages.

| Condition | 3 Months (young) | 18 Months (old) |
| --- | --- | --- |
| Naive control (`WT`) | `WT_5` | `WT_6` |
| *Sftpc*<sup>I73T</sup>, day 28 | `I73T_Week_4_5`, `I73T_Week_4_6` | `I73T_Week_4_18M_1`, `_2` |
| I73T + nerandomilast, day 28 | `Nerandomilast_3M` | `Nerandomilast_18M` |

Fibrosis is induced in the I73T model and assessed 4 weeks (day 28) later;
treated animals receive nerandomilast over the same window.

**Note on n:** the I73T groups have n = 2 per age; WT and nerandomilast groups
have n = 1. Treat between-condition differences as exploratory — this design
supports describing cell-state shifts, not formal per-sample statistics.

Sample-to-condition mapping lives in `sample_metadata.csv` (columns:
`file`, `type`, `sample`, `age`). Every script reads condition labels
from this file rather than parsing filenames — if you add a sample, edit the
manifest, not the code.

## Data availability

**The data is not in this repository, and never should be.** `.h5ad` and `.loom`
objects run to hundreds of MB each; Git stores every version of every file
forever, so committing them would permanently bloat the repo.

<!-- TODO: fill in once deposited -->
- Raw / processed data: GEO accession `GSE_______`
- Working copies live on the lab workstation under `lab/nerandomilast/`

To reproduce: download the samples into one folder, then set `DATA` in
`config.py` to point at it. That is the only path anything needs.

## Pipeline

Run in this order. Each step writes a file the next step reads.

| # | Step | Script | In → Out |
| --- | --- | --- | --- |
| 1 | Convert GEO looms | `loom_to_h5ad.py` | `*.loom` → `*.h5ad` |
| 2 | Filter + normalize | `preprocessing.ipynb` | per-sample `.h5ad` → `nerandomilast_preprocessed.h5ad` |
| 3 | Batch integration | `integration.py` | preprocessed → `nerandomilast_integrated.h5ad` |
| 4 | Compartments | `Annotation_code.ipynb` | integrated → `nerandomilast_annotated.h5ad` + one file per compartment |
| 5 | Fine annotation | *not yet in this repo* | each compartment → annotated with `cell_type` |
| 6 | Composition plots | `frequency_analysis.py` | annotated → stacked bar chart |
| 7 | Signalling trajectory | `cellchat_trajectory_plot.R` | CellChat strengths CSV → trajectory PDF/PNG |

**Steps 2, 4 and 5 are notebooks on purpose.** Each one requires looking at
output to decide what happens next — QC violins in step 2, marker genes and
clustering resolution in steps 4 and 5. Steps 1, 3, 6 and 7 are scripts because
they are pure transformations with no human in the loop.

### 1. Loom → h5ad (`loom_to_h5ad.py`)

GEO deposits arrive as a mix of `.loom` (velocyto) and `.h5ad`. This converts
looms to h5ad — faster to load, and one format downstream. Spliced, unspliced
and ambiguous layers are preserved, so RNA velocity remains possible later.

Each output is verified by reloading it and checking that shape, layer sums and
non-zero counts match the source exactly.

> ⚠️ **This script deletes the source `.loom` files and rewrites the manifest
> after a successful conversion.** Keep a backup of the raw downloads.

### 2. Preprocessing (`preprocessing.ipynb`)

Per sample, then concatenated:

- `min_genes = 200`, `min_cells = 3`
- gene counts trimmed to the **2nd–99th percentile** (per sample, so it adapts
  to differing sequencing depth rather than imposing one global cutoff)
- `pct_counts_mt < 10`
- `pct_counts_ribo < 20`
- raw counts preserved in `.layers['counts']` **before** normalization —
  scVI (step 3) and any DE test need integer counts, not log-normalized values
- `normalize_total(target_sum=1e4)` → `log1p`

Thresholds are the standard ones for mouse single-cell data.

Ends with QC violins per sample. **Look at them.** If one sample looks unlike
the others after filtering, that is worth knowing before it propagates.

### 3. Integration (`integration.py`)

3000 highly variable genes (`seurat_v3`, computed per sample so no single
sample dominates), then scVI with `batch_key='sample'`, default architecture
and training schedule. The latent representation is written to
`.obsm['X_scVI']` on the **full** gene set — HVG selection is used for training
only, so all genes stay available for marker inspection.

A GPU is recommended. CPU training works but is slow.

### 4. Compartments (`Annotation_code.ipynb`)

Neighbours on `X_scVI` → UMAP → Leiden at `resolution = 0.1` for broad
compartments (epithelial / mesenchymal / immune / endothelial), confirmed
against canonical markers:

| Marker | Compartment |
| --- | --- |
| `Epcam` | epithelial |
| `Pdgfra`, `Pdgfrb` | mesenchymal |
| `Ptprc` | immune |
| `Pecam1` | endothelial |
| `Msln` | mesothelial |
| `Mki67` | proliferating |

Each compartment is then written out as its own object for step 5.

Cluster numbers are tied to one specific scVI run — rerunning integration
renumbers them, and the cluster-to-compartment mapping has to be redone.

### 5. Fine annotation

Each compartment is re-clustered on its own at a higher resolution, marker
genes are inspected, and clusters are assigned a `cell_type` label.

**Not yet in this repo** — one notebook per compartment still to be added.

### 6–7. Downstream

`frequency_analysis.py` produces stacked bar charts of cell-type proportion per
condition. `cellchat_trajectory_plot.R` plots incoming vs outgoing signalling
strength for one cell type across conditions, drawing an arrow along
`WT → I73T_Week_4 → Nerandomilast` to show whether treatment returns signalling
toward baseline.

## Environment

```bash
# Python
pip install scanpy anndata leidenalg scvi-tools scikit-misc \
            h5py openpyxl pandas numpy scipy matplotlib
```

```r
# R
install.packages(c("ggplot2", "ggrepel"))
# CellChat: remotes::install_github("jinworks/CellChat")
```

<!-- TODO: pin versions once the analysis is final — scanpy and scvi-tools
     both change defaults between releases, which silently changes results.
     `pip freeze > requirements.txt` -->

## Repository layout

```
.
├── config.py                    # EDIT THIS: where the data lives
├── loom_to_h5ad.py              # 1. format conversion
├── preprocessing.ipynb          # 2. filter + normalize
├── integration.py               # 3. scVI
├── Annotation_code.ipynb        # 4. compartments
├── frequency_analysis.py        # 6. composition figures
├── cellchat_trajectory_plot.R   # 7. signalling trajectory
├── sample_metadata.csv          # sample → condition mapping
└── .gitignore                   # keeps data objects out of Git
```

## Known gaps

Honest list of what is not yet in this repo:

- **Fine annotation notebooks are missing** — step 5. The `cell_type` labels
  that every downstream figure depends on are not yet documented here.
- **The CellChat run itself is not here** — only the plotting script. The code
  that builds the CellChat objects lives on the lab workstation.
- **`loom_to_h5ad.py`, `frequency_analysis.py` and `cellchat_trajectory_plot.R`
  still have hard-coded paths** and have not been updated to use `config.py`.
- **`type` pools ages for WT but not for other conditions** (`WT` covers both
  3M and 18M, while `I73T_Week_4` and `I73T_Week_4_18M` are separate). Any
  `groupby('type')` therefore merges young and old controls. Group by
  `['type', 'age']` for age-aware comparisons.
- **Versions are unpinned.**
