# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress C Elegans (Cao et al 2017)
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
    fix_annotations,
    get_celltype_order,
    collect_gene_annotations,
    store_compressed_atlas,
    )


root_repo_folder = pathlib.Path(__file__).parent.parent.parent
atlas_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'c_elegans'
anno_fn = root_repo_folder / 'data' / 'gene_annotations' / 'c_elegans.PRJNA13758.WS287.annotations.gff3.gz'
webapp_fdn = root_repo_folder / 'webapp'
output_fdn = webapp_fdn / 'static' / 'atlas_data'
fn_out = output_fdn / 'c_elegans.h5'

species = 'c_elegans'

rename_dict = {
    'cell_types': {
        'Unclassified glia': 'glia',
        'gABAergic neurons': 'GABAergic neuron',
        'am/PH sheath cells': 'am/PH sheath',
        'pharyngeal epithelia': 'pharingeal epi',
        'distal tip cells': 'distal tip',
        'socket cells': 'socket',
        'excretory cells': 'excretory',
        'somatic gonad precursors': 'somatic gonad',
        'unclassified neurons': 'neuron',
        'dopaminergic neurons': 'dopaminergic neuron',
        'cholinergic neurons': 'cholinergic neuron',
        'ciliated sensory neurons': 'ciliated sensory neuron',
        'canal associated neurons': 'canal associated neuron',
        'touch receptor neurons': 'touch receptor neuron',
        'pharyngeal neurons': 'pharyngeal neuron',
        'oxygen sensory neurons': 'oxygen sensory neuron',
        'flp-1(+) interneurons': 'flp-1(+) interneuron',
        'other interneurons': 'neuron',
        'vulval precursors': 'vulval precursor',
        'coelomocytes': 'coelomocyte',
        'seam cells': 'seam',
        'sex myoblasts': 'sex myoblast',
    },
}

coarse_cell_types = [
]


celltype_order = [
    ('immune', [
        'glia',
    ]),
    ('epithelial', [
        'seam',
        'non-seam hypodermis',
        'pharyngeal epi',
        'coelomocyte',
        'distal tip',
    ]),
    ('endothelial', [
    ]),
    ('mesenchymal', [
        'body wall muscle',
        'pharyngeal muscle',
        'intestinal/rectal muscle',
        'sex myoblast',
        'am/PH sheath',
        'socket',
        'pharyngeal gland',
        'excretory',
        'rectum',
    ]),
    ('other', [
        'dopaminergic neuron',
        'cholinergic neuron',
        'ciliated sensory neuron',
        'GABAergic neuron',
        'canal associated neuron',
        'touch receptor neuron',
        'pharyngeal neuron',
        'oxygen sensory neuron',
        'flp-1(+) interneuron',
        'neuron',
        'germline',
        'somatic gonad',
        'vulval precursor',
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

    tissue_sources = get_tissue_data_dict(species, atlas_data_folder)
    tissues = list(tissue_sources.keys())
    for it, (tissue, full_atlas_fn) in enumerate(tissue_sources.items()):
        print(tissue)

        adata_tissue = anndata.read(full_atlas_fn)

        # Ignore cells with NaN in the cell.type column
        idx = adata_tissue.obs['cell.type'].isin(
                adata_tissue.obs['cell.type'].value_counts().index)
        adata_tissue = adata_tissue[idx]

        # cptt throughout
        sc.pp.normalize_total(
            adata_tissue,
            target_sum=1e4,
            key_added='coverage',
        )

        # Fix cell type annotations
        adata_tissue.obs['cell.type'] = adata_tissue.obs['cell.type'].str.lower()
        adata_tissue.obs['cellType'] = fix_annotations(
            adata_tissue, 'cell.type', 'c_elegans', tissue,
            rename_dict, coarse_cell_types,
        )

        # Age
        adata_tissue.obs['age'] = 'L2'

        # Correction might declare some cells as untyped/low quality
        # they have an empty string instead of an actual annotation
        if (adata_tissue.obs['cellType'] == '').sum() > 0:
            idx = adata_tissue.obs['cellType'] != ''
            adata_tissue = adata_tissue[idx]

        celltypes = get_celltype_order(
            adata_tissue.obs['cellType'].value_counts().index,
            celltype_order,
        )

        print('Add data to celltype group')
        genes = adata_tissue.var_names
        avg_ge = pd.DataFrame(
                np.zeros((len(genes), len(celltypes)), np.float32),
                index=genes,
                columns=celltypes,
                )
        frac_ge = pd.DataFrame(
                np.zeros((len(genes), len(celltypes)), np.float32),
                index=genes,
                columns=celltypes,
                )
        ncells_ge = pd.Series(
                np.zeros(len(celltypes), np.int64), index=celltypes,
                )
        for celltype in celltypes:
            idx = adata_tissue.obs['cellType'] == celltype
            Xidx = adata_tissue[idx].X
            avg_ge[celltype] = np.asarray(Xidx.mean(axis=0))[0]
            frac_ge[celltype] = np.asarray((Xidx > 0).mean(axis=0))[0]
            ncells_ge[celltype] = idx.sum()

        print('Add data to celltype-timepoint group')
        ages = adata_tissue.obs['age'].value_counts().index.tolist()
        ages.sort()
        columns_age = []
        for ct in celltypes:
            for age in ages:
                columns_age.append('_'.join([ct, 'Celegans', str(age)]))

        # Averages
        genes = adata_tissue.var_names
        avg_ge_tp = pd.DataFrame(
                np.zeros((len(genes), len(celltypes) * len(ages)), np.float32),
                index=genes,
                columns=columns_age,
                )
        frac_ge_tp = pd.DataFrame(
                np.zeros((len(genes), len(celltypes) * len(ages)), np.float32),
                index=genes,
                columns=columns_age,
                )
        ncells_ge_tp = pd.Series(
                np.zeros(len(columns_age), np.int64), index=columns_age,
                )
        for celltype in celltypes:
            adata_ct = adata_tissue[adata_tissue.obs['cellType'] == celltype]
            for age in ages:
                idx_age = (adata_ct.obs['age'] == age).values.nonzero()[0]
                if len(idx_age) == 0:
                    continue
                Xct_age = adata_ct.X[idx_age]
                label = '_'.join([celltype, 'TMC', str(age)])
                avg_ge_tp[label] = np.asarray(Xct_age.mean(axis=0))[0]
                frac_ge_tp[label] = np.asarray((Xct_age > 0).mean(axis=0))[0]
                ncells_ge_tp[label] = len(idx_age)

        compressed_atlas[tissue] = {
            'features': genes,
            'celltype': {
                'avg': avg_ge,
                'frac': frac_ge,
                'ncells': ncells_ge,
            },
            'celltype_dataset_timepoint': {
                'avg': avg_ge_tp,
                'frac': frac_ge_tp,
                'ncells': ncells_ge_tp,
            },
        }

    print('Consolidate gene list across tissues')
    needs_union = False
    genes = None
    for tissue, tdict in compressed_atlas.items():
        genest = list(tdict['features'])
        if genes is None:
            genes = genest
            continue
        if genest != genes:
            needs_union = True
            genes = set(genes) | set(genest)

    if needs_union:
        raise NotImplementedError('TODO: make union of features')

    print('Get gene annotations')
    gene_annos = collect_gene_annotations(anno_fn, genes)

    print('Store compressed atlas to file')
    store_compressed_atlas(
            fn_out,
            compressed_atlas,
            tissues,
            gene_annos,
            celltype_order,
    )
