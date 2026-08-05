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
| 11 | Composition plots | *not yet added* | annotated → stacked bar chart |
| 12 | Signalling trajectory | *not yet added* | CellChat strengths CSV → trajectory PDF/PNG |

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

Both combinations live in one script and run in one command:

```bash
python combine.py
```

They are in the same script because they are not independent — the second reads
`immune_annotated.h5ad`, which the first writes. The `COMBINATIONS` list at the
top of the file sets which objects go into which output; edit that if the file
names change.

It prints the cell and gene count of each input, then of the result, so the
arithmetic is visible: 16,778 + 5,587 + 12,387 = 34,752 immune cells, and
9,973 + 56,644 + 34,752 = 101,369 in the final object.

`join='outer'` keeps the union of genes. Each object was filtered on its own
after the split, so their gene sets differ slightly; an inner join would
quietly drop every gene missing from any one of them.

The final object carries 47 `cell_type` labels. No label is shared between
compartments, so the merge needs no relabelling — though `pDC2` does appear in
both the lymphocyte and the macrophage object, and those merge into one label.

### 11–12. Downstream

*Not yet added.* `frequency_analysis.py` produces stacked bar charts of
cell-type proportion per condition. `cellchat_trajectory_plot.R` plots incoming
vs outgoing signalling strength for one cell type across conditions, drawing an
arrow along `WT → I73T_Week_4 → Nerandomilast` to show whether treatment returns
signalling toward baseline.

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
├── .gitignore                   # keeps data objects out of Git
└── .gitattributes               # normalizes line endings across machines
```

Still to be added: `loom_to_h5ad.py` (1) and the downstream figure scripts
(11–12).

## Known gaps

Honest list of what is not yet in this repo:

- **The figure scripts are missing** — steps 11–12. The pipeline runs as far as
  `nerandomilast_annotated.h5ad`; nothing that reads it is here yet.
- **The CellChat run itself is not here** — only the plotting script is
  planned. The code that builds the CellChat objects lives on the lab
  workstation.
- **The scripts still to be added have hard-coded paths** and will need
  updating to use `config.py` before they go in.
- **`type` pools ages for WT but not for other conditions** (`WT` covers both
  3M and 18M, while `I73T_Week_4` and `I73T_Week_4_18M` are separate). Any
  `groupby('type')` therefore merges young and old controls. Group by
  `['type', 'age']` for age-aware comparisons.
- **Versions are unpinned.**
