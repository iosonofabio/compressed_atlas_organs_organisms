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
    fn_atlas = '../data/tabula_sapiens/TS_Heart.h5ad'
    adata = anndata.read_h5ad(fn_atlas)

    # NOTE: the human data is in some weird normalization between 0 and 10,
    # use adata.raw.X for computations to avoid log trasformations and whatnot

    print('Identify cell type column')
    # Find cell type column
    column = 'cell_ontology_class'
    # NOTE: there is a free annotation column that appears to be exactly the
    # same information yet again, and very much not "free"... ask Bob?

    print('Focus on macrophages for Toshie')
    amacro = adata[adata.obs[column].isin(['macrophage', 'dendritic cell'])]

    print('Back to raw counts')
    adata.X = adata.layers['raw_counts']

    if False:
        print('Convert to singlet')
        ds = singlet.Dataset.from_AnnData(amacro)
        
        print('Feature selection')
        features = ds.feature_selection.overdispersed_within_groups('donor', inplace=False)
        dsf = ds.query_features_by_name(features)

        # Manual features to consider
        # mouse: Siglecf, Fabp1, Cidec, Krt79, Ear6, Mcoln3, Rufy4, Acaa1b, Ear1, Il1rn, Vstm2a
        # human: SIGLEC14,FABP1,CIDEC,KRT79,RNASE2,MCOLN3,RUFY4,ACAA1,IL1RN,VSTM2A
        # NOTE: not all seem to work well in the heatmap 
        print('PCA')
        dsc = dsf.dimensionality.pca(n_dims=30, robust=False, return_dataset='samples')

        print('UMAP')
        vsu = dsc.dimensionality.umap()

        print('knn graph')
        edges = dsc.graph.knn(axis='samples', n_neighbors=10, return_kind='edges')

        print('Leiden clustering')
        ds.samplesheet['cluster'] = ds.cluster.leiden(
                axis='samples',
                edges=edges,
                resolution_parameter=0.0005,
                )
        nclusters = ds.samplesheet['cluster'].max() + 1

        print('Find markers')
        dsa = ds.average('samples', 'cluster')

        print('Plot dotplot')
        genes = [
            'MKI67',
            'CD52',
            'TIMP3',
            'IFI27',
            'CD163',
            'CDK1',
            'TREM2',
            'MARCO',
            'MAFF',
            'CTSW',
            'AQP3',
            'CD1C',
            'THBD',
            'CD68',
            'CD48',
            'CCL3',
            'CHIT1',
            'HSPA1A',
            'HSPA1B',
            'CA2',
            'RNASE6',
            'CA4',
            'C1QA',
            'C1QC',
            'PLAC8',
            'PLAUR',
            'CXCL2',
            'IL1RN',
            'KRT79',
            'FABP1',
            'CIDEC',
            'MCOLN3',
            'LPL',
            'LIPA',
            'ABCG1',
            'SIGLEC1',
            'IL18',
            'PSAP',
            'CSF1R',
            'CXCL2',
            'CXCL16',
            'CX3CR1',
            'IFITM3',
            'MS4A7',
            'PFN1',
            'GPX1',
            'GRN',
            'TYROBP',
            'IFNG',
            #'IFITM6',
            'IFITM2',
            #'MS4A6C',
            #'n_counts_UMIs',
            #'n_genes',
            'SFTPC',
            'PHLDA1',
            'HBEGF',
            'SNHG16',
            'CXCL3',

        ]
        fig, ax = plt.subplots()
        ds.plot.dot_plot(
            group_by='cluster',
            group_order=list(range(nclusters)),
            plot_list=genes,
            ax=ax,
            layout='vertical',
        )
        fig.tight_layout()

    else:
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
        sc.pl.umap(amacro, color=['CA4', 'AQP1', 'CD163', 'ID1', 'IFI27'])
        sc.pl.umap(amacro, color=['CD1C'])
        sc.pl.umap(amacro, color=column)
        sc.pl.umap(amacro, color='leiden')
        sc.pl.umap(amacro, color='donor')

        print('Find marker genes')
        sc.tl.rank_genes_groups(amacro, 'leiden', method='t-test')
        sc.pl.rank_genes_groups(amacro, n_genes=25, sharey=False)

        print('Limit to clear clusters')
        clus = ['0', '4', '9']
        amacro2 = amacro[amacro.obs['leiden'].isin(clus)]
        sc.tl.rank_genes_groups(amacro2, 'leiden', method='t-test')
        sc.pl.rank_genes_groups(amacro2, n_genes=25, sharey=False)

        print('Limit to DCs')
        clus = ['12', '13', '16']
        amacro2 = amacro[amacro.obs['leiden'].isin(clus)]
        sc.tl.rank_genes_groups(amacro2, 'leiden', method='t-test')
        sc.pl.rank_genes_groups(amacro2, n_genes=25, sharey=False)


        print('Plot dotplot')
        genes = [
            'TIMP3',
            'CAV1',
            'SPARC',
            'HIF1A',
            'EPAS1',
            'MKI67',
            'CD52',
            'TIMP3',
            'IFI27',
            'CD163',
            'CDK1',
            'TREM2',
            'MARCO',
            'MAFF',
            'CTSW',
            'AQP3',
            'CD1C',
            'THBD',
            'CD68',
            'CD48',
            'CCL3',
            'CHIT1',
            'HSPA1A',
            'HSPA1B',
            'CA2',
            'RNASE6',
            'CA4',
            'C1QA',
            'C1QC',
            'PLAC8',
            'PLAUR',
            'CXCL2',
            'IL1RN',
            'KRT79',
            'FABP1',
            'CIDEC',
            'MCOLN3',
            'LPL',
            'LIPA',
            'ABCG1',
            'SIGLEC1',
            'IL18',
            'PSAP',
            'CSF1R',
            'CXCL2',
            'CXCL16',
            'CX3CR1',
            'IFITM3',
            'MS4A7',
            'PFN1',
            'GPX1',
            'GRN',
            'TYROBP',
            'IFNG',
            #'IFITM6',
            'IFITM2',
            #'MS4A6C',
            #'n_counts_UMIs',
            #'n_genes',
            'SFTPC',
            'PHLDA1',
            'HBEGF',
            'SNHG16',
            'CXCL3',
        ]

        genes = [
            'CD1C',
            'MKI67',
            'CDK1',
            'ITGAE',
            'XCR1',
            'CD209',
            'MREG',
            'IL12B',
            'CD70',
            'CEACAM21',
        ]
        sc.pl.dotplot(amacro, genes, groupby='leiden', dendrogram=True)

        annodic = {
            '16': 'Proliferating DC',
            '12': 'DC I',
            '15': 'DC III',
            '13': 'DC II',
            '9': 'Alveolar mac',
            ('0', '4'): 'Interstitial mac',
        }
        annotations = amacro.obs[column].copy()
        annotations = annotations.cat.add_categories(list(annodic.values()))
        for key, val in annodic.items():
            if isinstance(key, str):
                key = [key]
            annotations.loc[amacro.obs['leiden'].isin(key)] = val

        #annotations.to_csv(
        #        '../data/tabula_sapiens/reannotations_macrophages_DCs.tsv',
        #        sep='\t',
        #        index=True,
        #        )

        plt.ion(); plt.show()


