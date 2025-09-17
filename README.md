# SPHERE

**Spatial Proximity-Weighted High-resolution Expression Reconstruction**  
A generalized bin-to-cell mapping framework for high-definition spatial transcriptomics (HD-ST)

![SPHERE Workflow](Figure1_workflow.pdf)

---

## 📖 Introduction
SPHERE is a Python package for bin-to-cell mapping in high-resolution spatial transcriptomics data (e.g., 10x Genomics Visium HD).  
It implements three mapping strategies:
- **Naïve** (full containment)
- **Weighted Area** (weighted by overlapping area)
- **SPHERE** (overlapping area × spatial proximity weighting)

SPHERE can directly convert bin-level expression matrices and cell segmentation results into a cell × gene expression matrix, enabling downstream clustering and biological analysis.

---

## 🔧 Installation

It is recommended to create a separate Conda environment to avoid dependency conflicts:

```bash
conda create -n SPHERE python=3.10 -y
conda activate SPHERE
```
Clone this repository:
```bash
git clone https://github.com/xiaoxiaoxier/SPHERE.git
cd SPHERE
```
Install dependencies:
```bash
pip install -r requirements.txt
```

Install SPHERE:
```bash
pip install -e .   # Development mode installation
```
## 📂 Data Preparation
SPHERE requires three types of inputs:

- bin_gdf: bin-aggregated spatial gene expression data (with geometry)
- cell_gdf: cell segmentation polygon data (with cell_id and geometry)
- ground truth (optional): real cell × gene expression matrix for evaluation
For a detailed data preparation workflow, see:
examples/synthetic_data_Xenium_HumanCRC20250916.ipynb
This notebook demonstrates how to generate the required files from Xenium public datasets.

## 🚀 Quick Start
```bash
from sphere import sphere_assignment, run_clustering
import pandas as pd, geopandas as gpd
from shapely import wkt

# Example patch ID
patch_id = "patch_0_0.csv"
base_dir = "/path/to/chunks"

bin_path = f"{base_dir}/bins_gdf/{patch_id}"
cell_path = f"{base_dir}/cells_gdf/{patch_id}"
gt_path   = f"{base_dir}/ground_truth_nuclei/{patch_id}"

# Read bin data
bin_df = pd.read_csv(bin_path)
bin_df["geometry"] = bin_df["geometry"].apply(wkt.loads)
bin_gdf = gpd.GeoDataFrame(bin_df, geometry="geometry")

# Read cell segmentation data
cell_df = pd.read_csv(cell_path)
cell_df["geometry"] = cell_df["geometry"].apply(wkt.loads)
cell_gdf = gpd.GeoDataFrame(cell_df, geometry="geometry")

# Read ground truth
gt = pd.read_csv(gt_path).set_index("cell_id")
gt = gt.apply(pd.to_numeric, errors='coerce').fillna(0)

# Find common genes
exclude_cols = {'assigned_bin_id', 'row', 'column', 'geometry', 'cell_id'}
common_genes = sorted(list(set(bin_gdf.columns) & set(gt.columns) - exclude_cols))

# Run SPHERE
pred_sphere = sphere_assignment(bin_gdf, cell_gdf, common_genes)

# Clustering
adata_clustered = run_clustering(pred_sphere)

```



## 📊 Downstream Analysis
run_clustering returns an AnnData object, which can be further analyzed using Scanpy:

Visualization: UMAP plots of clustering results
Marker gene analysis
Spatial visualization (with cell coordinates)
Example:
```bash
import scanpy as sc
sc.pl.umap(adata_clustered, color='leiden')
sc.tl.rank_genes_groups(adata_clustered, groupby='leiden', method='t-test')
sc.pl.rank_genes_groups(adata_clustered, n_genes=10)
```
## 📚 Examples
Full data preparation pipeline: examples/human_crc_prepare_inputs.ipynb
Single patch test: examples/test_patch_clustering.ipynb

## 📄 Citation
If you use SPHERE in your research, please cite:
```
@article{wang2025sphere,
  title={SPHERE: Spatial Proximity-Weighted High-resolution Expression Reconstruction for Visium High-Definition Transcriptomics},
  author={Wang, Zhuo and Li, Shangfu and Zhang, Chiping and Huang, Hsien-Da},
  journal={Bio},
  year={2025}
}
```
## 📜 License
This project is licensed under the MIT License.
