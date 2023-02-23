# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress Tabula Sapiens, heart.
'''
import os
import sys
import pathlib
import gzip
import h5py
import numpy as np
import pandas as pd

import anndata
import scanpy as sc

from utils import (
    get_tissue_data_dict,
    subannotate,
    correct_celltypes,
    get_celltype_order,
    )


root_repo_folder = pathlib.Path(__file__).parent.parent.parent
ts_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'tabula_sapiens'
anno_fn = root_repo_folder / 'data' / 'gene_annotations' / 'Homo_sapiens.GRCh38.109.gtf.gz'
webapp_fdn = root_repo_folder / 'webapp'
output_fdn = webapp_fdn / 'static' / 'atlas_data'
fn_out = output_fdn / 'tabula_sapiens.h5'


rename_dict = {
    'tissues': {
        'Large_Intestine': 'Colon',
    },
    'cell_types': {
    },
}

coarse_cell_types = [
]


celltype_order = [
    ('immune', [
    ]),
    ('epithelial', [
    ]),
    ('endothelial', [
    ]),
    ('mesenchymal', [
    ]),
    ('other', [
    ]),
    ('unknown', [
        'unknown',
    ])
]


if __name__ == '__main__':

    # Remove existing compressed atlas file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    compressed_atlas = {}

    tissue_sources = get_tissue_data_dict(ts_data_folder, rename_dict)
    tissues = list(tissue_sources.keys())
    for it, (tissue, full_atlas_fn) in enumerate(tissue_sources.items()):
        print(tissue)

        adata_tissue = anndata.read(full_atlas_fn)

        # Restart from raw data and renormalize
        adata_tissue = adata_tissue.raw.to_adata()

        # cptt throughout
        sc.pp.normalize_total(
            adata_tissue,
            target_sum=1e4,
            key_added='coverage',
        )

        # Fix cell type annotations
        adata_tissue.obs['cellType'] = correct_celltypes(
            adata_tissue, 'cell_ontology_class', 'human',
        )

        # Correction might declare some cells as untyped/low quality
        # they have an empty string instead of an actual annotation
        if (adata_tissue.obs['cellType'] == '').sum() > 0:
            idx = adata_tissue.obs['cellType'] != ''
            adata_tissue = adata_tissue[idx]

        celltypes = get_celltype_order(
            adata_tissue.obs['cellType'].value_counts().index,
            celltype_order,
        )

        break


    sys.exit()


    # Load data
    print('Load single cell data')
    fn_atlas = '../data/tabula_sapiens/TS_Heart.h5ad'
    adata = anndata.read_h5ad(fn_atlas)
    adata.uns['data_source'] = ts_source_url

    sys.exit()

    # We reannotated a few things ourselves
    print('Reannotate some cell types')
    correct_annotations(adata)

    # TODO: same for fibroblasts, etc.

    # NOTE: the human data is in some weird normalization between 0 and 10,
    # use adata.raw.X for computations to avoid log trasformations and whatnot

    print('Compress')
    # Find cell type column
    column = 'cell_ontology_class'
    # NOTE: there is a free annotation column that appears to be exactly the
    # same information yet again, and very much not "free"... ask Bob?

    # Set cell type order
    #celltypes = adata.obs[column].unique()
    celltypes = [
        'Mesothelial',
        'Adventitial FB',
        'Alveolar FB',
        'MyoF',
        'ASM',
        'VSM',
        'Pericyte',
        'Venous',
        'gCap',
        'Aerocyte',
        'Arterial II',
        'Lymphatic',
        'Bronchial EC',
        'Plasmablast',
        'B cell',
        'NK cell',
        'pDC',
        'DC I',
        'DC II',
        'DC III',
        'Alveolar mac',
        'Interstitial mac',
        'Monocyte',
        'Neutrophil',
        'Basophil',
        'Alveolar type I',
        'Alveolar type II',
        'Basal',
        'Ciliated',
        'Club',
        'Goblet',
        'Ionocyte',
        'Serous',
        'Mucous',
    ]

    # Merge some annotations together
    annotation_merges = {
        'Adventitial FB': ['adventitial cell'],
        'Alveolar FB': ['alveolar fibroblast'],
        'Mesothelial': ['mesothelial cell'],
        'MyoF': ['myofibroblast cell'],
        'ASM': ['bronchial smooth muscle cell'],
        'VSM': ['vascular associated smooth muscle cell'],
        'Pericyte': ['pericyte cell'],
        'Venous': ['vein endothelial cell'],
        'gCap': [
            'capillary endothelial cell',
            'lung microvascular endothelial cell',
        ],
        'Aerocyte': ['capillary aerocyte'],
        'Arterial II': ['endothelial cell of artery'],
        'Lymphatic': ['endothelial cell of lymphatic vessel'],
        'Bronchial EC': ['bronchial vessel endothelial cell'],
        'Ionocyte': ['pulmonary ionocyte'],
        'Plasmablast': ['plasma cell'],
        'NK cell': ['nk cell'],
        'IL cell': ['innate lymphoid cell'],
        'Bcell': ['b cell'],
        'T cell': [
            'T cell',
            'cd8-positive, alpha-beta t cell',
            'cd4-positive, alpha-beta t cell',
            'cd8-positive alpha-beta t cell',
            'cd4-positive alpha-beta t cell',
        ],
        'Alveolar type I': ['type i pneumocyte'],
        'Alveolar type II': ['type ii pneumocyte'],
        'Goblet': ['respiratory goblet cell'],
        'Ciliated': ['lung ciliated cell'],
        'Club': ['club cell'],
        'Basal': ['basal cell'],
        'Serous': ['serous cell of epithelium of bronchus'],
        'Mucous': ['respiratory mucous cell'],
        'DC I': ['DC I'],
        'DC II': ['DC II'],
        'DC III': ['DC III'],
        'Proliferating DC': ['Proliferating DC'],
        'pDC': ['plasmacytoid dendritic cell'],
        'Alveolar mac': ['Alveolar mac'],
        'Interstitial mac': ['Interstitial mac'],
        'Monocyte': [
            'classical monocyte',
            'non-classical monocyte',
            'intermediate monocyte',
            ],
        'Neutrophil': ['neutrophil'],
        'Basophil': ['basophil'],
    }
    adata.obs[column+'_compressed'] = ''
    for ct, csts in annotation_merges.items():
        adata.obs.loc[adata.obs[column].isin(csts), column+'_compressed'] = ct

    # Average, proportion expressing, number of cells
    genes = adata.var_names
    avg_exp = pd.DataFrame(
            np.zeros((len(genes), len(celltypes)), np.float32),
            index=genes,
            columns=celltypes)
    frac_exp = avg_exp.copy()
    ncells = pd.Series(np.zeros(len(celltypes), np.int64), index=celltypes)
    for ct in celltypes:
        # Avg gene expression (use adata.raw.X, compute cptt after averaging)
        avg_exp[ct] = np.asarray(adata[adata.obs[column+'_compressed'] == ct].raw.X.mean(axis=0))[0]
        avg_exp[ct] *= 1e4 / avg_exp[ct].sum()
        # Proportion expressing (use adata.raw.X)
        frac_exp[ct] = np.asarray((adata[adata.obs[column+'_compressed'] == ct].raw.X > 0).mean(axis=0))[0]
        # Number of cells
        ncells[ct] = (adata.obs[column+'_compressed'] == ct).sum()

    # Store to file
    print('Store compressed atlas')
    fn_out = output_fdn / 'human_compressed_heart_atlas.h5'
    with pd.HDFStore(fn_out) as h5_data:
        h5_data.put('gene_expression/celltype/average', avg_exp.T)
        h5_data.put('gene_expression/celltype/fraction', frac_exp.T)
        h5_data.put('gene_expression/celltype/cell_count', ncells)

    print('Compute gene friends')
    # Compute top correlated genes
    counts = avg_exp.values.T
    genes = avg_exp.index

    # focus on highly expressed genes
    idx = counts.max(axis=0) > 30
    counts = counts[:, idx]
    genes = genes[idx]

    ngenes = len(genes)
    mu = counts.mean(axis=0)
    countsc = counts - mu
    sigma = counts.std(axis=0)

    friends = {}
    for i in range(ngenes):
        if ((i+1) % 100) == 0:
            print(f'Gene {i+1} out of {ngenes}', end='\r')
        corr_i = (countsc[:, i] * counts.T).mean(axis=1) / (sigma[i] * sigma.T)
        corr_i[np.isnan(corr_i)] = 0
        # Self is always the highest
        itop = np.argsort(corr_i)[::-1][1:6]

        # Cutoff at 20%
        tmp = corr_i[itop]
        tmp = tmp[tmp >= 0.2]
        itop = itop[:len(tmp)]

        friends_i = genes[itop]
        friends[genes[i]] = np.array(friends_i)
    print(f'Gene {i+1} out of {ngenes}')

    print('Store gene friends')
    output_fn = output_fdn / 'human_gene_friends.h5'
    with h5py.File(output_fn, 'w') as f:
        h5g = f.create_group('gene_expression')
        for gene, friends_i in friends.items():
            lmax = max(len(x) for x in friends_i)
            friends_i = friends_i.astype('S'+str(lmax))
            dset = h5g.create_dataset(gene, data=friends_i)


