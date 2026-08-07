"""Concatenate annotated objects into one.

Runs each combination in order - the later ones read what the first writes:

  1. the three annotated immune lineages  -> immune_annotated.h5ad
  2. the three annotated compartments     -> nerandomilast_annotated.h5ad
  3. mesenchyme + immune                  -> mes_immune.h5ad

All paths are relative to DATA in config.py.

    python combine.py
"""
import anndata as ad
import scanpy as sc

from config import DATA

# order matters: immune_annotated.h5ad is written by the first combination and
# read by the two after it
COMBINATIONS = [
    (['lymphocytes_annotated.h5ad',
      'granulocytes_annotated.h5ad',
      'macrophages_annotated.h5ad'], 'immune_annotated.h5ad'),

    (['epithelial_annotated.h5ad',
      'mesenchyme_annotated.h5ad',
      'immune_annotated.h5ad'], 'nerandomilast_annotated.h5ad'),

    (['mesenchyme_annotated.h5ad',
      'immune_annotated.h5ad'], 'mes_immune.h5ad'),
]

for inputs, output in COMBINATIONS:
    print(f'--- {output} ---')

    adatas = []
    for name in inputs:
        a = sc.read_h5ad(DATA / name)
        print(f'{name}: {a.n_obs:,} cells, {a.n_vars:,} genes')
        adatas.append(a)

    # join='outer' keeps the union of genes. Each object was filtered on its
    # own after the split, so their gene sets differ slightly; an inner join
    # would quietly drop every gene missing from any one of them.
    adata = ad.concat(adatas, join='outer')

    print(f'combined: {adata.n_obs:,} cells, {adata.n_vars:,} genes')
    print(f"{adata.obs['cell_type'].nunique()} cell types")

    adata.write_h5ad(DATA / output)
    print(f'wrote {output}\n')
