# vim: fdm=indent
'''
author:     Fabio Zanini
date:       21/02/23
content:    Reannotate TMS pancreas
'''
import os
import sys
import pathlib

import numpy as np
import pandas as pd

import anndata
import scanpy as sc


if __name__ == '__main__':

    # Locate the h5ad file (this might differ for you)
    pancreas_h5ad_file = 'tabula-muris-senis-droplet-processed-official-annotations-Pancreas.h5ad'

    # Read file into RAM
    adata = anndata.read(pancreas_h5ad_file)

    # Restart from raw data and renormalize
    adata = adata.raw.to_adata()

    # cptt throughout, like Emily's RNA data
    sc.pp.normalize_total(
        adata,
        target_sum=1e4,
        key_added='coverage',
    )

    # log
    sc.pp.log1p(adata)

    # HVG
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    # PCA
    sc.pp.pca(adata)

    # graph
    sc.pp.neighbors(adata)

    # cluster the graph
    sc.tl.leiden(adata)

    # markers
    sc.tl.rank_genes_groups(adata, 'leiden')
    top_markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(30)

    # plot matrix of markers and clusters
    sc.pl.rank_genes_groups_matrixplot(adata)

    # and now... what the heck are those clusters?
    # are there any endothelial cells?
    # useful, non-pancreas-specific markers:
    # - Acta2/Mustn1/Mgp/Aspn: (smooth) muscle
    # - Cox4i2/Higd1b: pericytes
    # - Col6a2: mesenchymal
    # - Cdh5/Pecam1: endothelial
    # - Gja5/Bmx: arterial
    # - Vwf/Slc6a2: venous
    # - Prox1: lymphatics
    # - Cd34: usually immune
