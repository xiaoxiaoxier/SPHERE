# SPHERE
SPHERE: A generalized bin-to-cell mapping method for high-definition spatial transcriptomics


This repository provides:
- Preprocessing using [bin2cell](https://github.com/JEFworks/bin2cell) for cell segmentation
- Three bin-to-cell mapping methods: Naïve, Weighted Area, SPHERE (ours)
- Downstream analysis with Scanpy

## Installation
```bash
git clone https://github.com/<yourname>/SPHERE.git
cd SPHERE
pip install -r requirements.txt

## Quick Start
Example: Mouse Brain VisiumHD


from sphere.preprocessing import run_bin2cell_preprocessing
from sphere.assignment import sphere_assignment
from sphere.clustering import run_clustering

# Step 1: Preprocessing & segmentation
adata = run_bin2cell_preprocessing(
    path="path/to/binned_outputs",
    source_image_path="path/to/tissue_image.tif",
    spaceranger_image_path="path/to/spatial"
)

# Step 2: Mapping bins to cells using SPHERE
# (bin_gdf, cell_gdf, common_genes should come from your processed data)
pred_sphere = sphere_assignment(bin_gdf, cell_gdf, common_genes)

# Step 3: Clustering
adata_clustered = run_clustering(pred_sphere)


Examples
See 
examples/mouse_brain_demo.ipynb
 for a full walkthrough.

License
MIT License

Citation
If you use SPHERE in your research, please cite:

Your Name et al. SPHERE: A generalized bin-to-cell mapping approach for HD-ST, Journal, Year.


---

