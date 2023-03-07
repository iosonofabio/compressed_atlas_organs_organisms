# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress Tabula Sapiens.
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
        'cd24 neutrophil': 'neutrophil',
        'cd4-positive, alpha-beta t cell': 'T',
        'cd8-positive, alpha-beta t cell': 'T',
        'erythroid progenitor': 'erythroid',
        'nk cell': 'NK',
        'hematopoietic stem cell': 'HSC',
        'nampt neutrophil': 'neutrophil',
        'memory b cell': 'B',
        'naive b cell': 'B',
        'myeloid progenitor': 'myeloid',
        'plasmablast': 'plasma cell',
        'enterocyte of epithelium of large intestine': 'enterocyte',
        'immature enterocyte': 'enterocyte',
        'paneth cell of epithelium of large intestine': 'paneth',
        'mature enterocyte': 'enterocyte',
        'b cell': 'B',
        'large intestine goblet cell': 'goblet',
        'transit amplifying cell of large intestine': 'transit amp',
        'goblet cell': 'goblet',
        'intestinal crypt stem cell': 'crypt',
        'intestinal crypt stem cell of large intestine': 'crypt',
        'intestinal enteroendocrine cell': 'enteroendocrine',
        'gut endothelial cell': 'endothelial',
        'mast cell': 'mast',
        'intestinal tuft cell': 'brush',
        'cardiac muscle cell': 'cardiomyocyte',
        'cardiac endothelial cell': 'coronary',
        'fibroblast of cardiac tissue': 'fibroblast',
        'smooth muscle cell': 'smooth muscle',
        'cd4-positive helper t cell': 'T',
        'kidney epithelial cell': 'epithelial',
        'endothelial cell': 'endothelial',
        'type i pneumocyte': 'AT1 epi',
        'type ii pneumocyte': 'AT2 epi',
        'basal cell': 'basal',
        'classical monocyte': 'monocyte',
        'club cell': 'club',
        'non-classical monocyte': 'monocyte',
        'capillary endothelial cell': 'capillary',
        'respiratory goblet cell': 'goblet',
        'lung ciliated cell': 'ciliated',
        'capillary aerocyte': 'CAP2',
        'vein endothelial cell': 'venous',
        'lung microvascular endothelial cell': 'capillary',
        'adventitial cell': 'fibroblast',
        'dendritic cell': 'dendritic',
        'intermediate monocyte': 'monocyte',
        'pericyte cell': 'pericyte',
        'endothelial cell of artery': 'arterial',
        'cd4-positive alpha-beta t cell': 'T',
        'bronchial smooth muscle cell': 'smooth muscle',
        'vascular associated smooth muscle cell': 'vascular smooth muscle',
        'cd8-positive alpha-beta t cell': 'T',
        'endothelial cell of lymphatic vessel': 'lymphatic',
        'bronchial vessel endothelial cell': 'capillary',
        'pulmonary ionocyte': 'ionocyte',
        'plasmacytoid dendritic cell': 'plasmacytoid',
        'mesothelial cell': 'mesothelial',
        'serous cell of epithelium of bronchus': 'serous',
        'myofibroblast cell': 'smooth muscle',
        'respiratory mucous cell': 'mucous',
        'pancreatic acinar cell': 'acinar',
        'pancreatic ductal cell': 'ductal',
        'myeloid cell': 'myeloid',
        't cell': 'T',
        'pancreatic stellate cell': 'stellate',
        'pancreatic beta cell': 'beta',
        'pancreatic pp cell': 'PP',
        'pancreatic alpha cell': 'alpha',
        'pancreatic delta cell': 'delta',
        'epithelial cell': 'epithelial',
        'tongue muscle cell': 'striated muscle',
        'schwann cell': 'schwann',
    },
}

coarse_cell_types = [
    'endothelial',
    'immune cell',
]


celltype_order = [
    ('immune', [
        'HSC',
        'neutrophil',
        'basophil',
        'granulocyte',
        'mast',
        'myeloid',
        'monocyte',
        'macrophage',
        'dendritic',
        'erythroid',
        'erythrocyte',
        'B',
        'plasma cell',
        'T',
        'NK',
        'plasmacytoid',
    ]),
    ('epithelial', [
        'epithelial',
        'goblet',
        'brush',
        'crypt',
        'transit amp',
        'enterocyte',
        'paneth',
        'AT1 epi',
        'AT2 epi',
        'club',
        'ciliated',
        'ductal',
        'acinar',
        'keratinocyte',
        'basal',
        'serous',
        'mucous',
    ]),
    ('endothelial', [
        'arterial',
        'venous',
        'coronary',
        'capillary',
        'CAP2',
        'lymphatic',
    ]),
    ('mesenchymal', [
        'fibroblast',
        'alveolar fibroblast',
        'cardiomyocyte',
        'stellate',
        'striated muscle',
        'smooth muscle',
        'vascular smooth muscle',
        'pericyte',
        'mesothelial',
    ]),
    ('other', [
        'enteroendocrine',
        'hepatocyte',
        'ionocyte',
        'alpha',
        'beta',
        'PP',
        'delta',
        'schwann',
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
        adata_tissue.obs['cellType'] = fix_annotations(
            adata_tissue, 'cell_ontology_class', 'human', tissue,
            rename_dict, coarse_cell_types,
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
        # NOTE: see supplementary table 1 of the Science paper
        adata_tissue.obs['age'] = adata_tissue.obs['donor'].map({
            'TSP7': 69, 'TSP14': 59, 'TSP4': 38,
        })
        ages = adata_tissue.obs['age'].value_counts().index.tolist()
        ages.sort()
        columns_age = []
        for ct in celltypes:
            for age in ages:
                columns_age.append('_'.join([ct, 'TS', str(age)]))

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
                label = '_'.join([celltype, 'TS', str(age)])
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

    print('Add gene annotations')
    gene_annos = collect_gene_annotations(anno_fn, genes)

    print('Store compressed atlas to file')
    store_compressed_atlas(
            fn_out,
            compressed_atlas,
            tissues,
            gene_annos,
            celltype_order,
    )
