"""Stacked bar chart of cell type composition per condition.

Cross-tabulates cell_type against type, normalized within each condition, so
the bars show what a condition is made of rather than how many cells it has -
sample sizes differ, and raw counts would only show that.

Works on any annotated object. OBJECT, TITLE and COLORS all describe the same
object and are changed together:

    OBJECT   which .h5ad to read, relative to DATA in config.py
    TITLE    what goes above the chart
    COLORS   one entry per cell_type in that object

The script prints the cell types it found before plotting, so pointing OBJECT
somewhere new and running it once gives you the list to fill COLORS with.

    python frequency_analysis.py
"""
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt

from config import DATA

OBJECT = 'macrophages_annotated.h5ad'
TITLE = 'Macrophages Cell Type Frequency by Type'

# The set below belongs to the macrophage object. Point OBJECT at a different
# one and this whole dict is replaced by that object's cell types - a label
# with no entry raises KeyError where the column colours are looked up.
COLORS = {
    'Alveolar Macrophages': '#d62728',
    'Spp1+ Alveolar Macrophages (MoAMs)': '#7a7a7a',
    'Interstitial Macrophages': '#e9900b',
    'Spp1+ Interstitial Macrophages': '#f519a1',
    'Classical Monocytes': '#2ca02c',
    'Activated Monocytes': '#c6d128',
    'Cycling Myeloid Cells': '#21a7e6',
    'cDC1': '#4cb3c0',
    'cDC2': '#396b4a',
    'pDC1': '#8b4ad4',
    'pDC2': '#000000',
}

adata = sc.read_h5ad(DATA / OBJECT)

print(f'{OBJECT}: {adata.n_obs:,} cells')
print(f"{adata.obs['cell_type'].nunique()} cell types to colour:")
for label in sorted(adata.obs['cell_type'].unique()):
    print(f'    {label}')

# normalize='index' makes each row sum to 1, so every condition becomes a
# composition rather than a count
df = pd.crosstab(adata.obs['type'], adata.obs['cell_type'], normalize='index')
df_percent = (df * 100).round(1)

fig, ax = plt.subplots(figsize=(14, 8))
df_percent.plot(kind='bar', stacked=True, ax=ax,
                color=[COLORS[column] for column in df_percent.columns])

ax.legend(title='Cell type', bbox_to_anchor=(1.05, 1), loc='upper left')

ax.set_title(TITLE)
ax.set_xlabel('')
ax.set_ylabel('Proportion (%)')
ax.set_xticks(range(len(df_percent.index)))
ax.set_xticklabels(df_percent.index, rotation=45)
ax.grid(False)

fig.tight_layout()
plt.show()
