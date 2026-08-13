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
| 1 | Convert GEO looms | `loom_to_h5ad.py` *(not yet added)* | `*.loom` → `*.h5ad` |
| 2 | Filter + normalize | `preprocessing.ipynb` | per-sample `.h5ad` → `nerandomilast_preprocessed.h5ad` |
| 3 | Batch integration | `integration.py` | preprocessed → `nerandomilast_integrated.h5ad` |
| 4 | Compartments | `compartments.ipynb` | integrated → `nerandomilast_compartments.h5ad`, `epithelial.h5ad`, `mesenchyme.h5ad`, `immune.h5ad` |
| 5 | Epithelial annotation | `epithelial.ipynb` | `epithelial.h5ad` → `epithelial_annotated.h5ad` |
| 6 | Mesenchyme annotation | `mesenchyme.ipynb` | `mesenchyme.h5ad` → `mesenchyme_annotated.h5ad` |
| 7 | Immune sub-split | `immune_split.ipynb` | `immune.h5ad` → `lymphocytes.h5ad`, `granulocytes.h5ad`, `macrophages.h5ad` |
| 8 | Immune annotation | `lymphocytes.ipynb`, `granulocytes.ipynb`, `macrophages.ipynb` | each lineage → `*_annotated.h5ad` |
| 9 | Combine immune | `combine.py` | the three lineages → `immune_annotated.h5ad` |
| 10 | Combine compartments | `combine.py` | the three compartments → `nerandomilast_annotated.h5ad` |
| 11 | Composition plots | `frequency_analysis.py` | any annotated object → stacked bar chart |
| 12 | Subset for CellChat | `subset_cellchat_input.py` | annotated → `nerandomilast_10celltypes.h5ad` |
| 13 | Run CellChat | `cellchat_run.R` | subset → `cellchat_<condition>.rds` |
| 14 | Pathway figures | `cellchat_figures.R` | CellChat objects → circle + heatmap PDF per pathway set |

**The annotation steps are notebooks on purpose.** Each requires looking at
output to decide what happens next — QC violins in step 2, markers and
clustering resolution in steps 4–8. Steps 1, 3, 9–12 are scripts because they
are pure transformations with no human in the loop.

The immune compartment takes three steps rather than one because it is split
again — into macrophage, lymphocyte and granulocyte objects — each cleaned and
annotated on its own before being recombined. That sub-grouping is recorded in
`obs['immune_compartment']`.

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

### 4. Compartments (`compartments.ipynb`)

Neighbours on `X_scVI` → UMAP → Leiden at `resolution = 0.1`, deliberately
coarse. Each cluster is assigned to a compartment by reading canonical markers:

| Marker | Compartment |
| --- | --- |
| `Epcam` | Epithelial |
| `Pdgfra`, `Pdgfrb`, `Msln` | Mesenchyme |
| `Ptprc` | Immune |
| `Pecam1` | Endothelial |

Mesothelial cells (`Msln`) are kept inside the mesenchymal compartment rather
than split out. Clusters expressing none of these markers are labelled `other`.

**No cells are filtered in this step.** Three new objects are written —
epithelial, mesenchyme, immune — and the endothelial and `other` cells simply
are not among them:

- **Endothelial** (5,123 cells) is a scope decision. This project concerns the
  epithelial, mesenchymal and immune response to nerandomilast; the endothelium
  was not part of the question.
- **`other`** (234 cells) matched none of the four markers, so there is no
  compartment to carry them into.

The full object, all compartment labels intact, is kept as
`nerandomilast_compartments.h5ad` — the record of where every cell went,
including the ones that stop here.

Cluster numbers are tied to one specific scVI run — rerunning integration
renumbers them, and the cluster-to-compartment mapping has to be redone from the
markers.

### 5–8. Annotation

Every compartment follows the same two-pass shape:

1. **Cluster, inspect, remove junk.** Clusters are removed when they have no
   distinguishing marker of their own, or when they co-express markers of
   another compartment — residual `Pecam1` endothelium surviving the split is
   the usual offender, and this is where it gets caught. Mitochondrial and
   ribosomal filtering is not repeated here; that happened per sample in step 2.
2. **Re-cluster what remains, then annotate.** Removing cells changes every
   neighbourhood, so the graph and UMAP are rebuilt before the clustering the
   labels are assigned to.

**The clustering resolutions are not principled values.** They are tuned by eye
— raised until the cells of interest separate into a cluster of their own, and
stopped before real populations fragment. What works for one compartment says
nothing about the next.

Across the three compartments this removed a further 6,708 cells.

The process was iterative: compartments were re-clustered and inspected over
several passes, and cluster numbering changed between them. **The per-pass
cluster numbers are therefore not recoverable, and this repo does not claim to
regenerate them.** What is fixed is the outcome — the `cell_type` labels on the
final object, which is what every downstream figure reads.

**Step 5, epithelial (`epithelial.ipynb`)** — 9,973 cells, 9 cell types.
Alveolar identity from `Sftpc` (AT2), `Ager`/`Hopx` (AT1), with `Sftpc` +
`Lcn2` marking Activated AT2 and `Cldn4` marking Intermediate Cells. The
Transitional AT2-AT1 state is read from the combined marker picture rather than
any single gene. Airway identity from `Scgb1a1` (secretory), `Muc5b` (goblet),
`Foxj1` (ciliated), `Trp63` (basal).

**Step 6, mesenchyme (`mesenchyme.ipynb`)** — 56,644 cells, 10 cell types. Each
population is called from a panel rather than a single gene, because fibroblast
subsets differ by degree more than by presence or absence: adventitial,
alveolar, fibrotic, inflammatory and peribronchial fibroblasts, plus cycling
fibroblasts, pericytes, smooth muscle and mesothelial cells.

Follicular dendritic cells are the exception. They sit inside the adventitial
fibroblasts — the perivascular niche where tertiary lymphoid structures form —
and only separate once the resolution is raised. `Cxcl13` is the call, and
`Tnfsf13b` (BAFF) is read alongside it, the same axis this project follows into
the FDC–B cell signalling analysis.

**Step 7, immune sub-split (`immune_split.ipynb`)** — no annotation happens
here. `resolution = 0.1`, coarse enough to separate three megagroups and
nothing finer: lymphocytes (`Cd3e`, `Nkg7`, `Cd19`, `Xbp1`), granulocytes
(`S100a9`, `Csf1`) and macrophages (`Marco`, `Ccl2`, `Spp1`, `Cd74`, `Atox1`).
Each is written out so it can be re-clustered at a resolution suited to it.

**Step 8, immune annotation** — one notebook per lineage, each following the
same two-pass shape.

- `lymphocytes.ipynb` — 16,778 cells, 13 cell types. T subsets (`Cd3e`, `Cd4`,
  `Cd8a`, `Foxp3`, `Trdc`), B and plasma (`Cd19`, `Xbp1`), NK (`Nkg7`), plus
  ILC2s and the cycling fractions.
- `granulocytes.ipynb` — 5,587 cells, 5 cell types. Neutrophils (`Retnlg`,
  `S100a9`) and eosinophils (`Csf1`), each split into a pair of sibling
  clusters, plus mast cells.
- `macrophages.ipynb` — 12,387 cells, 11 cell types. Despite the name this
  object holds the whole myeloid lineage: resident and monocyte-derived
  macrophages, monocytes, and the dendritic cell subsets.

### 9–10. Combining (`combine.py`)

All the combinations live in one script and run in one command:

```bash
python combine.py
```

| Output | Built from |
| --- | --- |
| `immune_annotated.h5ad` | lymphocytes + granulocytes + macrophages |
| `nerandomilast_annotated.h5ad` | epithelial + mesenchyme + immune |
| `mes_immune.h5ad` | mesenchyme + immune |

They are in one script because they are not independent — the last two both
read `immune_annotated.h5ad`, which the first writes, so the order matters. The
`COMBINATIONS` list at the top of the file sets which objects go into which
output; edit that to add or rename one.

It prints the cell and gene count of each input, then of the result, so the
arithmetic is visible: 16,778 + 5,587 + 12,387 = 34,752 immune cells, and
9,973 + 56,644 + 34,752 = 101,369 in the final object.

`join='outer'` keeps the union of genes. Each object was filtered on its own
after the split, so their gene sets differ slightly; an inner join would
quietly drop every gene missing from any one of them.

The final object carries 47 `cell_type` labels. No label is shared between
compartments, so the merge needs no relabelling — though `pDC2` does appear in
both the lymphocyte and the macrophage object, and those merge into one label.

### 11. Composition plots (`frequency_analysis.py`)

Stacked bar chart of cell type composition per condition. The cross-tabulation
uses `normalize='index'`, so each bar sums to 100% — the chart shows what a
condition is *made of*, not how many cells it has. That matters here: sample
sizes differ several-fold, so raw counts would mostly show sequencing depth.

Runs on any annotated object. Three constants at the top describe the same
object and change together:

| Constant | |
| --- | --- |
| `OBJECT` | which `.h5ad` to read, relative to `DATA` |
| `TITLE` | what goes above the chart |
| `COLORS` | one entry per `cell_type` in that object |

Colours are assigned explicitly rather than left to matplotlib so a cell type
keeps the same colour across every figure it appears in. The trade-off is that
`COLORS` must cover every label — a missing one raises `KeyError`.

To point it at a different compartment, set `OBJECT` and run it once: it prints
the cell types it found before plotting, which is the list `COLORS` needs.

```bash
python frequency_analysis.py
```

### 12. Subset for CellChat (`subset_cellchat_input.py`)

CellChat scales badly with the number of cell types, and most of the 47 in the
annotated object have nothing to do with the question this analysis asks. This
cuts it down to the **FDC–B cell axis** — 50,410 cells, ten cell types:

| Stromal | B lineage |
| --- | --- |
| Adventitial Fibroblast | B Cells |
| Alveolar Fibroblast | Naive B Cells |
| Fibrotic Fibroblast | Plasma Cells |
| Inflammatory Fibroblast | Cycling Plasma Cells |
| FDC | pDC2 |

Which is what BAFF signalling needs: FDC are the source, B cells the target.

Dropped are every T, NK, myeloid and epithelial population, and the mesenchymal
types outside the fibroblast lineage — pericytes, smooth muscle, peribronchial
fibroblasts, mesothelial cells, cycling fibroblasts.

Nothing is re-clustered or re-labelled. `cell_type` and `leiden` carry through
unchanged, so a cell keeps the identity it was given during annotation.

```bash
python subset_cellchat_input.py
```

### 13. Run CellChat (`cellchat_run.R`)

Reads the subset through `zellkonverter` and builds one CellChat object per
condition. **Inference only — it produces no figures.** Plotting is a separate
step, so the expensive part runs once and any number of figures can be drawn
from the saved objects afterwards.

```bash
Rscript cellchat_run.R
```

**Resume-safe.** Any `cellchat_<condition>.rds` already on disk is skipped, and
the `sce.rds` conversion is cached. A run that dies partway can be restarted
without redoing the expensive parts — which matters, because
`computeCommunProb` is slow.

Two settings change the results and are worth knowing before comparing to
anything else:

| | |
| --- | --- |
| `MAX_CELLS_PER_TYPE = 800` | populations above this are randomly downsampled |
| `MIN_CELLS = 10` | a cell type below this **in a given condition** is dropped from that condition only |
| `NBOOT = 20` | permutations for the significance test; CellChat's default is 100 |

`MIN_CELLS` is why some conditions end up with nine cell types rather than ten:
FDC are rare, and in the sparser conditions they fall under the threshold.
`liftCellChat` pads both objects to their union before any pairwise comparison,
so the mismatch does not break the comparisons.

`set.seed(42)` makes the downsampling reproducible. `N_WORKERS` should be set
to the core count of whatever machine this runs on.

Output is `cellchat_<condition>.rds`, one per condition, plus the cached
`sce.rds`. Everything downstream reads those.

### 14. Pathway figures (`cellchat_figures.R`)

Circle and heatmap figures for each entry in `SETS` — either a single pathway
or several summed together:

```r
SETS <- list(
  TNF            = "TNF",
  RANKL          = "RANKL",
  TWEAK          = "TWEAK",
  BAFF           = "BAFF",
  LIGHT          = "LIGHT",
  TNFsuperfamily = c("TNF", "RANKL", "TWEAK", "BAFF", "LIGHT")
)
```

```bash
Rscript cellchat_figures.R
```

Two PDFs per set: `circle_<set>_5conditions.pdf` with the five conditions side
by side, and `heatmap_<set>_5conditions.pdf` with one panel per condition.

**Scales are shared within a set, not across sets.** One `edge.weight.max` and
one colour ramp apply to all five conditions, so thickness and intensity mean
the same thing in every panel — the point of putting them side by side. Sharing
across sets instead would make the weaker pathways invisible next to `TNF`.

Not every member of a set is detected in every condition. Each panel sums only
what that condition has, and its title records which — so
`TNFsuperfamily - WT (TNF+TWEAK+BAFF+LIGHT)` next to
`TNFsuperfamily - Nerandomilast (TNF+RANKL+TWEAK+BAFF)` is telling you the two
panels are not summing the same members. Conditions where nothing is detected
are drawn as a labelled blank so the five columns stay aligned across figures.

**`CELLCHAT_DIR` is set inside this script, not in `config.py`.** R cannot read
the Python config, so the path to the CellChat objects is defined separately.

## Cell counts

Where cells were lost between integration and the final annotated object:

| Stage | Cells | Change |
| --- | --- | --- |
| After preprocessing and integration | 113,434 | — |
| Endothelial not carried forward | 108,311 | −5,123 |
| `other` not carried forward | 108,077 | −234 |
| Junk, doublets and residual endothelium removed per compartment | **101,369** | −6,708 |

No cells are removed at the compartment split itself — endothelial and `other`
are simply not carried into the three objects that go forward. All quality
filtering after integration happens inside the individual compartment objects,
where each is re-clustered at a resolution fine enough to separate junk from
real populations.

12,065 cells in total, spread proportionally across all 8 samples (873–2,670
each, tracking sample size) — the losses are not concentrated in any one sample.

The 101,369 remaining cells carry 47 `cell_type` labels across three
compartments, and are the basis of every figure below.

## Environment

```bash
# Python
pip install scanpy anndata leidenalg scvi-tools scikit-misc \
            h5py pandas numpy scipy matplotlib
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

What is in the repo now:

```
.
├── README.md
├── config.py                    # EDIT THIS: where the data lives
├── sample_metadata.csv          # sample → condition mapping
├── preprocessing.ipynb          # 2. filter + normalize
├── integration.py               # 3. scVI
├── compartments.ipynb           # 4. compartments
├── epithelial.ipynb             # 5. epithelial annotation
├── mesenchyme.ipynb             # 6. mesenchyme annotation
├── immune_split.ipynb           # 7. immune sub-split
├── lymphocytes.ipynb            # 8. lymphocyte annotation
├── granulocytes.ipynb           # 8. granulocyte annotation
├── macrophages.ipynb            # 8. myeloid annotation
├── combine.py                   # 9, 10. concatenate annotated objects
├── frequency_analysis.py        # 11. composition figures
├── subset_cellchat_input.py     # 12. cut down to the FDC-B cell axis
├── cellchat_run.R               # 13. EDIT PATHS: build CellChat objects
├── cellchat_figures.R           # 14. EDIT CELLCHAT_DIR: circle + heatmap figures
└── .gitignore                   # keeps data objects out of Git
```

Still to be added: `loom_to_h5ad.py` (1).

## Known gaps

Honest list of what is not yet in this repo:

- **`loom_to_h5ad.py` is missing** — step 1. Everything from preprocessing
  onward is here; the format conversion that feeds it is not.
- **`type` pools ages for WT but not for other conditions** (`WT` covers both
  3M and 18M, while `I73T_Week_4` and `I73T_Week_4_18M` are separate). Any
  `groupby('type')` therefore merges young and old controls. Group by
  `['type', 'age']` for age-aware comparisons.
- **Versions are unpinned.**
