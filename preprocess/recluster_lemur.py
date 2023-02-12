# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress Tabula Sapiens.
'''
import os
import sys
import numpy as np
import pandas as pd
import h5py
import anndata

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet/build/lib')
import singlet
import scanpy as sc


data_fdn = '../webapp/static/scData/'


if __name__ == '__main__':

    # Load data
    print('Load single cell data')
    fn_atlas = '../data/tabula_microcebus/Lung_FIRM_hvg.h5ad'
    adata = anndata.read_h5ad(fn_atlas)

    # NOTE: the lemur data is in some weird normalization between 0 and 8.79,
    # and adata.raw is None. Hmm, there must be a logp1 there, let's try to
    # undo that transformation by hand. After np.expm1 the sum of each cell is
    # 10000, so there you go, multiply by 100 and done.
    adata.X = 100 * np.expm1(adata.X)

    print('Identify cell type column')
    # Find cell type column
    column = 'cell_ontology_class_v1'
    column2 = 'free_annotation_v1'

    annos = []

    print('Aerocytes are already annotated')
    aerocytes = adata.obs_names[(adata.obs[column2] == 'capillary aerocyte cell')]
    aerocytes = pd.Series(
        index=aerocytes,
        data=['Aerocyte' for x in aerocytes],
    )
    annos.append(aerocytes)

    print('Same with alveolar/adventitial fibros')
    # Actually it's ok to just use:
    # fibroblast of lung -> alveolar
    # fibroblast -> adventitial (good enough for now)

    print('Dendritic cells need reclustering? Not really')
    dcs = adata.obs.loc[adata.obs[column].isin([
        'conventional dendritic cell',
        'dendritic cell',
    ]), 'free_annotation_v1']
    dcs = (dcs.replace('conventional dendritic cell', 'DC I')
              .replace('mature dendritic cell', 'DC III')
              .replace('dendritic cell (FLT3+ IGSF6+)', 'DC II'))
    annos.append(dcs)

    if False:
        amacro = adata[adata.obs[column].isin([
            'conventional dendritic cell',
            'dendritic cell',
        ])]

        print('Normalize')
        sc.pp.normalize_total(amacro, target_sum=1e6)

        print('Get log')
        sc.pp.log1p(amacro)

        print('Feature selection')
        sc.pp.highly_variable_genes(amacro, min_mean=0.0125, max_mean=3, min_disp=0.5)
        amacro.raw = amacro
        amacro = amacro[:, amacro.var.highly_variable]

        print('PCA')
        sc.tl.pca(amacro, svd_solver='arpack')

        print('knn graph')
        sc.pp.neighbors(amacro, n_neighbors=10, n_pcs=40)

        print('UMAP embedding')
        sc.tl.umap(amacro)

        print('Cluster')
        sc.tl.leiden(amacro)

        print('Plot a few UMAPs')
        sc.pl.umap(amacro, color=['leiden', column, 'free_annotation_v1'])
        sc.pl.umap(amacro, color=['XCR1', 'MREG', 'IL12B'])
        sc.pl.umap(amacro, color=column)
        sc.pl.umap(amacro, color='donor')

        print('Find marker genes')
        sc.tl.rank_genes_groups(amacro, 'free_annotation_v1', method='t-test')
        sc.pl.rank_genes_groups(amacro, n_genes=25, sharey=False)


    annotations = pd.concat(annos)
    annotations.to_csv(
            '../data/tabula_microcebus/reannotations.tsv',
            sep='\t',
            index=True,
            )
