# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress Tabula Microcebus.
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
tmc_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'tabula_microcebus'
anno_fn = root_repo_folder / 'data' / 'gene_annotations' / 'Microcebus_murinus.Mmur_3.0.109.gtf.gz'
webapp_fdn = root_repo_folder / 'webapp'
output_fdn = webapp_fdn / 'static' / 'atlas_data'
fn_out = output_fdn / 'tabula_microcebus.h5'


rename_dict = {
    'tissues': {
    },
    'cell_types': {
        'erythroid lineage cell': 'erythroid',
        'unassigned': 'unknown',
        'erythroid progenitor cell': 'erythroid',
        'CD4-positive, alpha-beta T cell': 'T',
        'CD8-positive, alpha-beta T cell': 'T',
        'B cell': 'B',
        'hematopoietic precursor cell': 'hematopoietic',
        'natural killer cell': 'NK',
        'granulocyte monocyte progenitor cell': 'granulocytopoietic',
        'mature NK T cell': 'NK',
        'T cell': 'T',
        'myeloid cell': 'unknown',  # TODO
        'plasmacytoid dendritic cell': 'plasmacytoid',
        'megakaryocyte progenitor cell': 'megakaryocyte-erythroid',
        'dendritic cell': 'dendritic',
        'endothelial cell of sinusoid': 'capillary',
        'conventional dendritic cell': 'dendritic',
        'fat cell': 'adipocyte',
        'enterocyte of epithelium of large intestine': 'enterocyte',
        'epithelial cell of large intestine': 'epithelial',
        'mesothelial cell': 'mesothelial',
        'large intestine goblet cell': 'goblet',
        'intestinal enteroendocrine cell': 'enteroendocrine',
        'cell': 'unknown',
        'vascular associated smooth muscle cell': 'vascular smooth muscle',
        'endothelial cell': 'capillary',
        'regular ventricular cardiac myocyte': 'cardiomyocyte',
        'pericyte cell': 'pericyte',
        'regular atrial cardiac myocyte': 'cardiomyocyte',
        'endothelial cell of lymphatic vessel': 'lymphatic',
        'Purkinje myocyte': 'cardiomyocyte',
        'mesothelial cell of epicardium': 'mesothelial',
        'nodal myocyte': 'cardiomyocyte',
        'vasa recta ascending limb cell': 'capillary',
        'vasa recta descending limb cell': 'capillary',
        'kidney proximal convoluted tubule epithelial cell': 'proximal tubule epi',
        'kidney loop of Henle thin descending limb epithelial cell': 'Henle limb epi',
        'kidney loop of Henle thick ascending limb epithelial cell': 'Henle limb epi',
        'kidney loop of Henle thin ascending limb epithelial cell': 'Henle limb epi',
        'kidney proximal straight tubule epithelial cell': 'proximal tubule epi',
        'renal alpha-intercalated cell': 'intercalated',
        'renal beta-intercalated cell': 'intercalated',
        'glomerular endothelial cell': 'glomerular',
        'capillary endothelial cell': 'capillary',
        'epithelial cell of proximal tubule': 'proximal tubule epi',
        'renal principal cell': 'principal',
        'myofibroblast cell': 'myofibroblast',
        'kidney distal convoluted tubule epithelial cell': 'distal tubule epi',
        'macula densa epithelial cell': 'macula densa',
        'kidney collecting duct cell': 'collecting duct epi',
        'renal intercalated cell': 'intercalated',
        'innate lymphoid cell': 'innate lymphoid',
        'kidney loop of Henle epithelial cell': 'Henle limb epi',
        'reticular cell': 'reticular',
        'non-myelinating Schwann cell': 'schwann',
        'type II pneumocyte': 'AT2 epi',
        'fibroblast of lung': 'fibroblast',
        'type I pneumocyte': 'AT1 epi',
        'endothelial cell of artery': 'arterial',
        'vein endothelial cell': 'venous',
        'lung ciliated cell': 'ciliated',
        'brush cell of bronchus': 'brush',
        'club cell': 'club',
        'basal cell of epithelium of bronchus': 'basal',
        'myelinating Schwann cell': 'schwann',
        'pancreatic acinar cell': 'acinar',
        'pancreatic ductal cell': 'ductal',
        'pancreatic B cell': 'beta',
        'pancreatic A cell': 'alpha',
        'pancreatic D cell': 'delta',
        'pancreatic PP cell': 'PP',
        'epithelial cell of exocrine pancreas': 'epithelial',
        'oral mucosa squamous cell': 'squamous',
        'skeletal muscle satellite stem cell': 'satellite',
        'tendon cell': 'tendon',
        'fast muscle cell': 'striated muscle',
    },
}

coarse_cell_types = [
    'endothelial',
    'immune cell',
    'lymphocyte',
]


celltype_order = [
    ('immune', [
        'HSC',
        'hematopoietic',
        'neutrophil',
        'basophil',
        'eosinophil',
        'granulocytopoietic',
        'granulocyte',
        'mast',
        'myeloid',
        'monocyte',
        'alveolar macrophage',
        'macrophage',
        'dendritic',
        'megakaryocyte-erythroid',
        'erythroid',
        'erythrocyte',
        'platelet',
        'B',
        'plasma cell',
        'T',
        'NK',
        'plasmacytoid',
        'innate lymphoid',
    ]),
    ('epithelial', [
        'epithelial',
        'goblet',
        'brush',
        'crypt',
        'transit amp',
        'enterocyte',
        'paneth',
        'proximal tubule epi',
        'distal tubule epi',
        'podocyte',
        'Henle limb epi',
        'collecting duct epi',
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
        'squamous',
        'intercalated',
        'principal',
        'macula densa',
        'urothelial',
    ]),
    ('endothelial', [
        'arterial',
        'venous',
        'coronary',
        'capillary',
        'CAP2',
        'lymphatic',
        'glomerular',
    ]),
    ('mesenchymal', [
        'fibroblast',
        'alveolar fibroblast',
        'myofibroblast',
        'cardiomyocyte',
        'stellate',
        'tendon',
        'satellite',
        'striated muscle',
        'smooth muscle',
        'vascular smooth muscle',
        'pericyte',
        'mesothelial',
        'reticular',
        'preosteoblast',
        'osteoblast',
        'adipocyte',
    ]),
    ('other', [
        'neuron',
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

    tissue_sources = get_tissue_data_dict(tmc_data_folder, rename_dict)
    tissues = list(tissue_sources.keys())
    for it, (tissue, full_atlas_fn) in enumerate(tissue_sources.items()):
        print(tissue)

        adata_tissue = anndata.read(full_atlas_fn)

        # There is no raw data, but the actual data is log1p of cptt
        adata_tissue.X = np.expm1(adata_tissue.X)
        adata_tissue.obs['coverage'] = adata_tissue.obs['nCount_RNA']

        # Fix cell type annotations
        adata_tissue.obs['cellType'] = fix_annotations(
            adata_tissue, 'cell_ontology_class_v1', 'mouselemur', tissue,
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
