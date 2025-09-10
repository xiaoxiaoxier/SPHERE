# SPHERE
SPHERE: A generalized bin-to-cell mapping method for high-definition spatial transcriptomics
# SPHERE

**SPHERE**: Spatial Probabilistic Hybrid Expression Reconstruction for high-definition spatial transcriptomics.

This repository provides:
- Preprocessing using [bin2cell](https://github.com/JEFworks/bin2cell) for cell segmentation
- Three bin-to-cell mapping methods: Naïve, Weighted Area, SPHERE (ours)
- Downstream analysis with Scanpy

## Installation
```bash
git clone https://github.com/<yourname>/SPHERE.git
cd SPHERE
pip install -r requirements.txt
