# SPHERE

**Spatial Proximity-Weighted High-resolution Expression Reconstruction**  
A generalized bin-to-cell mapping framework for high-definition spatial transcriptomics (HD-ST)

![SPHERE Workflow](docs/SPHERE_workflow.svg)

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

git clone https://github.com/xiaoxiaoxier/SPHERE.git
cd SPHERE