"""
Batch integration with scVI.

Set DATA in config.py before running. GPU is recommended; training runs on
CPU but is very slow.
"""
import scanpy as sc
import scvi

from config import DATA

adata = sc.read_h5ad(DATA / 'nerandomilast_preprocessed.h5ad')
print(f'loaded {adata.shape}')

# 3000 HVGs, computed per sample so that no single sample dominates the
# selection. seurat_v3 expects raw integer counts, hence layer='counts'.
sc.pp.highly_variable_genes(
    adata, n_top_genes=3000, flavor='seurat_v3',
    batch_key='sample', layer='counts'
)

# Train on the HVG subset only, but write the latent representation back onto
# the full object so every gene stays available for marker analysis later.
n_before = adata.n_vars
adata_hvg = adata[:, adata.var['highly_variable']].copy()

scvi.model.SCVI.setup_anndata(adata_hvg, layer='counts', batch_key='sample')

model = scvi.model.SCVI(adata_hvg)
model.train()

adata.obsm['X_scVI'] = model.get_latent_representation()
assert adata.n_vars == n_before  # subsetting must not have altered the full object

out = DATA / 'nerandomilast_integrated.h5ad'
adata.write_h5ad(out)
print(f'wrote {out}')
