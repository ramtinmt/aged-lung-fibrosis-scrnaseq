"""Machine-specific paths.

This is the only file you edit when moving to a different computer.
Every script imports DATA from here, so the path is written down once.
"""
from pathlib import Path

# -----------------------------------------------------------------
# EDIT THIS: the folder holding the .h5ad files downloaded from GEO
# -----------------------------------------------------------------
DATA = Path("/path/to/nerandomilast/data")

# Repo root — where sample_metadata.csv lives. Leave this alone.
REPO = Path(__file__).resolve().parent
