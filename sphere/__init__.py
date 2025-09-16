"""
SPHERE package

Spatial Proximity-Weighted High-resolution Expression Reconstruction
for high-definition spatial transcriptomics (HD-ST).

This package provides:
- Bin-to-cell mapping methods (Naïve, Weighted Area, SPHERE)
- Cell × gene expression matrix clustering
- Preprocessing for bin2cell segmentation & data preparation

Usage example:
--------------
from sphere import sphere_assignment, run_clustering

pred_sphere = sphere_assignment(bin_gdf, cell_gdf, common_genes)
adata_clustered = run_clustering(pred_sphere)
"""

# 公共接口导入
from .assignment import naive_assignment, weighted_by_area, sphere_assignment
from .clustering import run_clustering
from .preprocessing import run_bin2cell_preprocessing

__all__ = [
    "naive_assignment",
    "weighted_by_area",
    "sphere_assignment",
    "run_clustering",
    "run_bin2cell_preprocessing",
]