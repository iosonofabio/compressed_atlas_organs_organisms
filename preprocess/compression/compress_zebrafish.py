# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress zebrafish (Wagmer et al 2018)
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
atlas_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'danio_rerio'
anno_fn = root_repo_folder / 'data' / 'gene_annotations' / 'Danio_rerio.GRCz11.109.gtf.gz'
webapp_fdn = root_repo_folder / 'webapp'
output_fdn = webapp_fdn / 'static' / 'atlas_data'
fn_out = output_fdn / 'd_rerio.h5'

species = 'd_rerio'

rename_dict = {
    'cell_types': {
        '24hpf-optic cup': 'optic cup',
        '24hpf-pharyngeal arch - cd248b': '',
        '24hpf-epidermal - col7a1l': '',
        '24hpf-pharyngeal pouch': '',
        '24hpf-pharyngeal arch - ndnf': 'ph arches',
        '24hpf-pharyngeal arch - tbx1': 'ph arches',
        '24hpf-pharyngeal arch - lbx1a': 'ph arches',
        '24hpf-hatching gland': 'hatching gland',
        '24hpf-periderm': 'peridermal',
        '24hpf-notocord': 'notocord',
        '24hpf-otic placode': 'otic placode',
        '24hpf-epidermal - olfactory placode': 'olfactory placode',
        '24hpf-lens': 'lens',
        '24hpf-erythroid': 'erythroid',
        '24hpf-macrophage': 'macrophage',
        '24hpf-leukocyte': 'leukocyte',
        '24hpf-pancreas primordium': 'pancreas primordium',
        '24hpf-neural crest - iridoblast': 'iridophore',
        '24hpf-neural crest - melanoblast': 'melanophore',
        '24hpf-neural crest - xanthophore': 'xantophore',
        '24hpf-pronephric duct': 'kidney epi',
        '24hpf-retina pigmented epithelium': 'retinal',
        '24hpf-pectoral fin bud': 'fin epi',
        '24hpf-epidermal - anterior': 'anterior epi',
        '24hpf-lateral line - krt15': 'epithelial',
        '24hpf-epidermal - rbp4': 'epithelial',
        '24hpf-epidermal - and1': 'epithelial',
        '24hpf-epidermal - kera': 'epithelial',
        '24hpf-epidermal - prr15la': 'epithelial',
        '24hpf-epidermal - atp1a1a.2': 'epithelial',
        '24hpf-epidermal - muc5ac': 'epithelial',
        '24hpf-epidermal - grhl3': 'epithelial',
        '24hpf-epidermal - acbd7': 'epithelial',
        '24hpf-epidermal - s100a11': 'epithelial',
        '24hpf-mesoderm - emp2': 'fibroblast',
        '24hpf-heart': 'cardiomyocyte',
        '24hpf-heart - hoxd9a': 'cardiomyocyte',
        '24hpf-heart - mature': 'cardiomyocyte',
        '24hpf-muscle - myl10': 'striated muscle',
        '24hpf-muscle - myl1': 'striated muscle',
        '24hpf-myotome': 'striated muscle',
        '24hpf-proctodeum': 'rectum',
        '24hpf-ionocyte - ca2': 'ionocyte',
        '24hpf-neural crest': 'neural crest',
        '24hpf-neural crest - mcamb': 'neural crest',
        '24hpf-neural crest - grem2': 'neural crest',
        '24hpf-neural - floorplate': 'neuron',
        '24hpf-neural - diencephalon posterior': 'neuron',
        '24hpf-differentiating neurons - sst1.1': 'neuron',
        '24hpf-neural - midbrain': 'neuron',
        '24hpf-neural - ventral hindbrain': 'neuron',
        '24hpf-neural - dorsal hindbrain': 'neuron',
        '24hpf-differentiating neurons': 'neuron',
        '24hpf-differentiating neurons - hmx': 'neuron',
        '24hpf-differentiating neurons - phox2a': 'neuron',
        '24hpf-neural - hindbrain roofplate': 'neuron',
        '24hpf-differentiating neurons - eomesa': 'neuron',
        '24hpf-differentiating neurons - dlx': 'neuron',
        '24hpf-neural - telencephalon': 'neuron',
        '24hpf-neural - hindbrain gsx1': 'neuron',
        '24hpf-neural - diencephalon ': 'neuron',
        '24hpf-neural - midbrain ventral nkx6.2': 'neuron',
        '24hpf-neural - posterior ventral nkx6.2': 'neuron',
        '24hpf-differentiating neurons - rohon beard': 'neuron',
        '24hpf-endoderm': 'endoderm',
        '24hpf-endothelial': 'capillary',
        '24hpf-endothelial - posterior': 'capillary',
        '24hpf-neural - dorsal spinal cord': 'spinal cord',
        '24hpf-tailbud - spinal cord': 'spinal cord',
        '24hpf-germline': 'germline',
        '24hpf-tailbud - PSM': 'PSM',
    },
}

coarse_cell_types = [
]


celltype_order = [
    ('immune', [
        'macrophage',
        'leukocyte',
        'erythroid',
    ]),
    ('epithelial', [
        'anterior epi',
        'fin epi',
        'ionocyte',
        'epithelial',
        'peridermal',
        'kidney epi',
        'retinal',
        'xantophore',
        'melanophore',
        'iridophore',
        'notocord',
    ]),
    ('endothelial', [
        'endoderm',
        'capillary',
    ]),
    ('mesenchymal', [
        'fibroblast',
        'striated muscle',
        'cardiomyocyte',
        'rectum',
        'pancreas primordium',
        'otic placode',
        'olfactory placode',
    ]),
    ('other', [
        'neuron',
        'spinal cord',
        'neural crest',
        'lens',
        'optic cup',
        'ph arches',
        'germline',
        'PSM',
        'hatching gland',
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
        idx = adata_tissue.obs['CellType'].isin(
                adata_tissue.obs['CellType'].value_counts().index)
        adata_tissue = adata_tissue[idx]

        # cptt throughout
        sc.pp.normalize_total(
            adata_tissue,
            target_sum=1e4,
            key_added='coverage',
        )

        # Fix cell type annotations
        adata_tissue.obs['cellType'] = fix_annotations(
            adata_tissue, 'CellType', 'c_elegans', tissue,
            rename_dict, coarse_cell_types,
        )

        # Age
        adata_tissue.obs['age'] = '24hpf'

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
