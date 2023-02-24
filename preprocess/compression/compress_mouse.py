# vim: fdm=indent
'''
author:     Carsten Knutsen, Fabio Zanini
date:       06/12/22
content:    Compress mouse data.

This comes from two sources:
- Emily's own data: this is compressed according to cell types and subtypes,
  which are annotated consistently across RNA- and ATAC-Seq.
- Tabula Muris Senis: this is only RNA-Seq, and does not use the exact same
  annotations as Emily's data. We use northstar and treasuremap to transfer the
  annotations from Emily's since they did a thorough job at annotating.
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
tms_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'tabula_muris_senis'
anno_fn = root_repo_folder / 'data' / 'gene_annotations' / 'mm10.ncbiRefSeq.gtf.gz'
webapp_fdn = root_repo_folder / 'webapp'
output_fdn = webapp_fdn / 'static' / 'atlas_data'
fn_out = output_fdn / 'tabula_muris_senis.h5'


rename_dict = {
    'tissues': {
        'Large_Intestine': 'Colon',
        'Heart_and_Aorta': 'Heart',
        'Marrow': 'Bone Marrow',
    },
    'cell_types': {
        'endothelial cell of coronary artery': 'coronary',
        'fibroblast of cardiac tissue': 'fibroblast',
        'endocardial cell': 'endocardial',
        'smooth muscle cell': 'smooth muscle',
        'cardiac neuron': 'neuron',
        'mast cell': 'myeloid',
        # See separate chat about vent VS atrial cardiomyocytes in TMS/Emily's
        'cardiomyocyte': 'ventricular',
        'precursor B cell': 'precursor B',
        'immature B cell': 'immature B',
        'late pro-B cell': 'late pro-B',
        'naive B cell': 'B',
        'naive T cell': 'T',
        'B cell': 'B',
        'T cell': 'T',
        'NK cell': 'NK',
        'enterocyte of epithelium of large intestine': 'enterocyte',
        'intestinal crypt stem cell': 'crypt',
        'epithelial cell of large intestine': 'epithelial',
        'large intestine goblet cell': 'goblet',
        'hematopoietic stem cell': 'HSC',
        'hematopoietic precursor cell': 'hematopoietic',
        'granulocytopoietic cell': 'granulocytopoietic',
        'megakaryocyte-erythroid progenitor cell': 'megakaryocyte-erythroid',
        'erythroid progenitor': 'erythroid',
        'kidney proximal convoluted tubule epithelial cell': 'proximal tubule epi',
        'epithelial cell of proximal tubule': 'proximal tubule epi',
        'kidney proximal straight tubule epithelial cell': 'proximal tubule epi',
        'kidney loop of Henle thick ascending limb epithelial cell': 'Henle limb epi',
        'kidney loop of Henle ascending limb epithelial cell': 'Henle limb epi',
        'kidney collecting duct principal cell': 'collecting duct epi',
        'kidney collecting duct epithelial cell': 'collecting duct epi',
        'kidney distal convoluted tubule epithelial cell': 'distal tubule epi',
        'brush cell': 'brush',
        'kidney cortex artery cell': 'arterial',
        'kidney mesangial cell': 'mesangial',
        'kidney capillary endothelial cell': 'capillary',
        'kidney cell': 'unknown',
        'fenestrated cell': 'fenestrated',
        'lung neuroendocrine cell': 'neuroendocrine',
        'classical monocyte': 'monocyte',
        'bronchial smooth muscle cell': 'smooth muscle',
        'intermediate monocyte': 'monocyte',
        'fibroblast of lung': 'alveolar fibroblast',
        'lung macrophage': 'macrophage',
        'non-classical monocyte': 'monocyte',
        'CD8-positive, alpha-beta T cell': 'T',
        'CD4-positive, alpha-beta T cell': 'T',
        'adventitial cell': 'unknown',
        'mature NK T cell': 'NKT',
        'vein endothelial cell': 'venous',
        'myeloid dendritic cell': 'dendritic',
        'pulmonary interstitial fibroblast': 'fibroblast',
        'type II pneumocyte': 'AT2 epi',
        'regulatory T cell': 'Treg',
        'smooth muscle cell of the pulmonary artery': 'vascular smooth muscle',
        'plasmacytoid dendritic cell': 'plasmacytoid',
        'pericyte cell': 'pericyte',
        'dendritic cell': 'dendritic',
        'endothelial cell of lymphatic vessel': 'lymphatic',
        'ciliated columnar cell of tracheobronchial tree': 'ciliated',
        'club cell of bronchiole': 'club',
        'pancreatic A cell': 'alpha',
        'pancreatic B cell': 'beta',
        'pancreatic D cell': 'delta',
        'pancreatic PP cell': 'PP',  # NOTE: ex-gamma cells
        'pancreatic stellate cell': 'stellate',
        'pancreatic acinar cell': 'acinar',
        'pancreatic ductal cel': 'ductal',  # NOTE: typo in the original
        'basal cell of epidermis': 'basal',
        'Langerhans cell': 'Langerhans',
        'leukocyte': 'macrophage',  # NOTE: kidney only
    }
}

coarse_cell_types = [
    'endothelial cell',  # pancreas
    'lymphocyte',  # kidney -> legit mixture with low-qual as well
    #'leukocyte',  # kidney
]

# We store an organism-wide complete ordering of cell types, and each tissue
# will cherry pick the necessary
celltype_order = [
    ('immune', [
        'HSC',
        'hematopoietic',
        'neutrophil',
        'basophil',
        'granulocytopoietic',
        'granulocyte',
        'promonocyte',
        'myeloid',
        'monocyte',
        'alveolar macrophage',
        'macrophage',
        'dendritic',
        'Langerhans',
        'megakaryocyte-erythroid',
        'proerythroblast',
        'erythroblast',
        'erythroid',
        'erythrocyte',
        'precursor B',
        'late pro-B',
        'immature B',
        'B',
        'plasma cell',
        'T',
        'Treg',
        'NKT',
        'NK',
        'plasmacytoid',
    ]),
    ('epithelial', [
        'epithelial',
        'goblet',
        'brush',
        'crypt',
        'enterocyte',
        'proximal tubule epi',
        'distal tubule epi',
        'podocyte',
        'Henle limb epi',
        'collecting duct epi',
        'AT2 epi',
        'club',
        'ciliated',
        'ductal',
        'acinar',
        'keratinocyte',
        'basal',
    ]),
    ('endothelial', [
        'arterial',
        'venous',
        'coronary',
        'fenestrated',
        'capillary',
        'lymphatic',
    ]),
    ('mesenchymal', [
        'fibroblast',
        'alveolar fibroblast',
        'endocardial',
        'ventricular',
        'stellate',
        'smooth muscle',
        'vascular smooth muscle',
        'pericyte',
        'mesangial',
    ]),
    ('other', [
        'neuron',
        'neuroendocrine',
        'alpha',
        'beta',
        'PP',
        'delta',
    ]),
    ('unknown', [
        'unknown',
    ]),
]


if __name__ == '__main__':

    # Remove existing compressed atlas file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    compressed_atlas = {}

    tissue_sources = get_tissue_data_dict(tms_data_folder, rename_dict)
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
            adata_tissue, 'cell_ontology_class', 'mouse',
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
        ages.sort(key=lambda x: int(x[:-1]))
        columns_age = []
        for ct in celltypes:
            for age in ages:
                columns_age.append('_'.join([ct, 'TMS', age]))

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
                label = '_'.join([celltype, 'TMS', age])
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
