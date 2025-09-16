"""
assignment.py
Bin-to-cell mapping methods for high-definition spatial transcriptomics.

Includes:
- naive_assignment
- weighted_by_area
- sphere_assignment
"""

import geopandas as gpd
import numpy as np


def naive_assignment(bin_gdf, cell_gdf, common_genes):
    """
    Assign bins to cells if fully contained (Naïve method).

    Parameters
    ----------
    bin_gdf : GeoDataFrame
        Bin-level polygons with gene expression columns.
    cell_gdf : GeoDataFrame
        Cell-segmentation polygons (must contain 'cell_id').
    common_genes : list
        List of genes to include in assignment.

    Returns
    -------
    DataFrame
        cell × gene matrix.
    """
    sj = gpd.sjoin(
        gpd.GeoDataFrame(bin_gdf, geometry="geometry"),
        gpd.GeoDataFrame(cell_gdf, geometry="geometry"),
        how="inner", predicate="within"
    )
    sj["cell_id"] = sj["index_right"].map(cell_gdf["cell_id"])
    return sj.groupby("cell_id")[common_genes].sum()


def weighted_by_area(bin_gdf, cell_gdf, common_genes):
    """
    Assign bins to cells proportional to geometric overlap area (Weighted Area method).
    """
    sj = gpd.sjoin(
        gpd.GeoDataFrame(bin_gdf, geometry="geometry"),
        gpd.GeoDataFrame(cell_gdf, geometry="geometry"),
        how="inner", predicate="intersects"
    )

    def safe_intersection(g1, g2):
        try:
            return g1.intersection(g2).area
        except Exception:
            return 0

    sj["overlap_area"] = sj.apply(
        lambda r: safe_intersection(r["geometry"],
                                    cell_gdf.loc[r["index_right"], "geometry"]),
        axis=1
    )
    sj["bin_area"] = sj["geometry"].area
    sj["weight"] = sj["overlap_area"] / sj["bin_area"]

    sj["cell_id"] = sj["index_right"].map(cell_gdf["cell_id"])
    weighted = sj[common_genes].multiply(sj["weight"], axis=0)
    weighted["cell_id"] = sj["cell_id"]
    return weighted.groupby("cell_id")[common_genes].sum()


def sphere_assignment(bin_gdf, cell_gdf, common_genes, alpha=0.01):
    """
    Assign bins to cells using SPHERE (overlap × proximity weighting).
    """
    sj = gpd.sjoin(
        gpd.GeoDataFrame(bin_gdf, geometry="geometry"),
        gpd.GeoDataFrame(cell_gdf, geometry="geometry"),
        how="inner", predicate="intersects"
    )

    def safe_intersection(g1, g2):
        try:
            return g1.intersection(g2).area
        except Exception:
            return 0

    sj["overlap_area"] = sj.apply(
        lambda r: safe_intersection(r["geometry"],
                                    cell_gdf.loc[r["index_right"], "geometry"]),
        axis=1
    )
    sj["bin_area"] = sj["geometry"].area
    sj["area_weight"] = sj["overlap_area"] / sj["bin_area"]

    # 距离权重
    sj["bin_cx"] = sj["geometry"].apply(lambda g: g.centroid.x)
    sj["bin_cy"] = sj["geometry"].apply(lambda g: g.centroid.y)
    cell_gdf["cx"] = cell_gdf["geometry"].apply(lambda g: g.centroid.x)
    cell_gdf["cy"] = cell_gdf["geometry"].apply(lambda g: g.centroid.y)
    sj["cell_cx"] = sj["index_right"].map(cell_gdf["cx"])
    sj["cell_cy"] = sj["index_right"].map(cell_gdf["cy"])
    sj["dist_weight"] = np.exp(-alpha * np.sqrt(
        (sj["bin_cx"] - sj["cell_cx"])**2 + (sj["bin_cy"] - sj["cell_cy"])**2
    ))

    # 最终权重
    sj["raw_w"] = sj["area_weight"] * sj["dist_weight"]
    sj["weight"] = sj.groupby(sj.index)["raw_w"].transform(lambda x: x / x.sum())

    # 按权重加和
    sj["cell_id"] = sj["index_right"].map(cell_gdf["cell_id"])
    weighted = sj[common_genes].multiply(sj["weight"], axis=0)
    weighted["cell_id"] = sj["cell_id"]
    return weighted.groupby("cell_id")[common_genes].sum()