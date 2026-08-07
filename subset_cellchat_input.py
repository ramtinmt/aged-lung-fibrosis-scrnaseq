"""Subset the annotated object to the cell types used in the CellChat analysis.

Keeps the populations relevant to this part of the project and drops the rest.
Nothing is re-clustered or re-labelled.

    python subset_cellchat_input.py
"""
import scanpy as sc

from config import DATA

SOURCE = 'nerandomilast_annotated.h5ad'
OUTPUT = 'nerandomilast_10celltypes.h5ad'

KEEP = [
    'Adventitial Fibroblast',
    'Alveolar Fibroblast',
    'Fibrotic Fibroblast',
    'Inflammatory Fibroblast',
    'FDC',
    'B Cells',
    'Naive B Cells',
    'Plasma Cells',
    'Cycling Plasma Cells',
    'pDC2',
]

adata = sc.read_h5ad(DATA / SOURCE)
print(f'{SOURCE}: {adata.n_obs:,} cells, {adata.obs["cell_type"].nunique()} cell types')

missing = [label for label in KEEP if label not in set(adata.obs['cell_type'])]
if missing:
    print(f'not present in {SOURCE}: {missing}')

adata = adata[adata.obs['cell_type'].isin(KEEP)].copy()
adata.obs['cell_type'] = adata.obs['cell_type'].cat.remove_unused_categories()

print(f'\n{adata.n_obs:,} cells kept, {adata.obs["cell_type"].nunique()} cell types')
print(adata.obs['cell_type'].value_counts().to_string())

adata.write_h5ad(DATA / OUTPUT)
print(f'\nwrote {OUTPUT}')
