import bin2cell as b2c
import scanpy as sc
import os

def run_bin2cell_preprocessing(path,
                               source_image_path,
                               spaceranger_image_path,
                               mpp=0.5,
                               outdir="bin2cell_output",
                               run_stardist=False):
    """
    Run basic bin2cell pre-processing to obtain cell segmentation masks.
    Parameters
    ----------
    path : str
        Path to the binned outputs from Visium HD.
    source_image_path : str
        Path to original tissue image (TIFF).
    spaceranger_image_path : str
        Path to spaceranger 'spatial' folder.
    mpp : float
        Microns per pixel for scaling.
    outdir : str
        Directory to save intermediate results.
    run_stardist : bool
        Whether to run stardist segmentation (else assumes npz labels exist).
    Returns
    -------
    AnnData
        AnnData object with segmentation labels.
    """
    os.makedirs(outdir, exist_ok=True)

    adata = b2c.read_visium(path,
                            source_image_path=source_image_path,
                            spaceranger_image_path=spaceranger_image_path)
    adata.var_names_make_unique()

    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_counts=1)

    b2c.scaled_he_image(adata, mpp=mpp, save_path=f"{outdir}/he.tiff")
    b2c.destripe(adata)

    if run_stardist:
        b2c.stardist(image_path=f"{outdir}/he.tiff",
                     labels_npz_path=f"{outdir}/he.npz",
                     stardist_model="2D_versatile_he",
                     prob_thresh=0.01)

    b2c.insert_labels(adata,
                      labels_npz_path=f"{outdir}/he.npz",
                      basis="spatial",
                      spatial_key="spatial_cropped_150_buffer",
                      mpp=mpp,
                      labels_key="labels_he")

    b2c.expand_labels(adata,
                      labels_key='labels_he',
                      expanded_labels_key="labels_he_expanded")

    return adata