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
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import anndata
import scanpy as sc


root_repo_folder = pathlib.Path(__file__).parent.parent.parent  
tms_data_folder = root_repo_folder / 'data' / 'full_atlases' / 'tabula_muris_senis'
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
        'endothelial cell of coronary artery': 'coronary cap',
        'fibroblast of cardiac tissue': 'fibroblast',
        'endocardial cell': 'endocardial',
        'smooth muscle cell': 'smooth muscle',
        'cardiac neuron': 'neuron',
        'mast cell': 'myeloid',
        # See separate chat about vent VS atrial cardiomyocytes in TMS/Emily's
        'cardiomyocyte': 'ventricular',
        'naive B cell': 'B cell',
        'naive T cell': 'T cell',
        'enterocyte of epithelium of large intestine': 'enterocyte',
        'intestinal crypt stem cell': 'crypt',
        'epithelial cell of large intestine': 'epithelial',
        'large intestine goblet cell': 'goblet',
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
    }
}

# We store an organism-wide complete ordering of cell types, and each tissue
# will cherry pick the necessary
celltype_order = [
    ('immune', [
        'hematopoietic',
        'basophil',
        'granulocytopoietic',
        'granulocyte',
        'promonocyte',
        'myeloid',
        'monocyte',
        'macrophage',
        'megakaryocyte-erythroid',
        'erythroid',
        'proerythroblast',
        'erythroblast',
        'precursor B cell',
        'late pro-B cell',
        'immature B cell',
        'B cell',
        'plasma cell',
        'T cell',
        'NK cell',
        'lymphocyte',  # TODO WTF kidney
        'leukocyte',  # TODO WTF kidney
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
    ]),
    ('endothelial', [
        'arterial',
        'coronary cap',
        'fenestrated cap',
        'capillary',
    ]),
    ('mesenchymal', [
        'fibroblast',
        'endocardial',
        'ventricular',
        'smooth muscle',
        'mesangial',
    ]),
    ('other', [
        'neuron',
    ]),
    ('unknown', [
        'unknown',
    ]),
]


def get_tissue_data_dict(atlas_folder):
    '''Get a dictionary with tissue order and files'''
    result = []

    fns = os.listdir(atlas_folder)
    fns = [x for x in fns if '.h5ad' in x]

    for fn in fns:
        tissue = fn.split('-')[-1].split('.')[0]
        tissue = rename_dict['tissues'].get(tissue, tissue)
        result.append({
            'tissue': tissue,
            'filename': atlas_folder / fn,
        })

    result = pd.DataFrame(result).set_index('tissue')['filename']

    # Order tissues alphabetically
    result = result.sort_index()

    return result


def correct_celltypes(celltypes_raw):
    '''Correct cell types in each tissue according to known dict'''
    celltypes_new = celltypes_raw.copy()
    for ctraw, celltype in rename_dict['cell_types'].items():
        celltypes_new[celltypes_new == ctraw] = celltype
    return celltypes_new


def get_celltype_order(celltypes_unordered):
    '''Use global order to reorder cell types for this tissue'''
    celltypes_ordered = []
    for broad_type, celltypes_broad_type in celltype_order:
        for celltype in celltypes_broad_type:
            if celltype in celltypes_unordered:
                celltypes_ordered.append(celltype)

    missing_celltypes = False
    for celltype in celltypes_unordered:
        if celltype not in celltypes_ordered:
            print('Missing celltype:', celltype)
            missing_celltypes = True

    if missing_celltypes:
        raise IndexError("Missing cell types!")

    return celltypes_ordered


if __name__ == '__main__':

    # Remove existing compressed atlas file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    tissue_sources = get_tissue_data_dict(tms_data_folder)
    for it, (tissue, full_atlas_fn) in enumerate(tissue_sources.items()):
        if it < 4:
            continue
        print(tissue)

        adata_tissue = anndata.read(full_atlas_fn)

        # Restart from raw data and renormalize
        adata_tissue = adata_tissue.raw.to_adata()

        # cptt throughout, like Emily's RNA data
        sc.pp.normalize_total(
            adata_tissue,
            target_sum=1e4,
            key_added='coverage',
        )

        adata_tissue.obs['cellType'] = correct_celltypes(
            np.asarray(adata_tissue.obs['cell_ontology_class'])
        )
        celltypes = get_celltype_order(
            adata_tissue.obs['cellType'].value_counts().index,
        )

        continue

    sys.exit()

    if False:
        print('Add data to celltype group')
        genes = adata_tissue.var_names
        avg_exp = pd.DataFrame(
                np.zeros((len(genes), len(celltypes)), np.float32),
                index=genes,
                columns=celltypes,
                )
        frac_exp = pd.DataFrame(
                np.zeros((len(genes), len(celltypes)), np.float32),
                index=genes,
                columns=celltypes,
                )
        ncells_exp = pd.Series(
                np.zeros(len(celltypes), np.int64), index=celltypes,
                )

        #  TODO

        print('Add data to celltype-timepoint group')
        ages_tms = adata_tissue.obs['age'].value_counts().index.tolist()
        ages_tms.sort(key=lambda x: int(x[:-1]))
        columns_age = []
        for ct in celltypes:
            for age in ages_tms:
                columns_age.append('_'.join([ct, 'TMS', age]))

        # Averages
        genes = adata_tissue.var_names
        avg_exp_tp_tms = pd.DataFrame(
                np.zeros((len(genes), len(celltypes) * len(ages_tms)), np.float32),
                index=genes,
                columns=columns_age,
                )
        frac_exp_tp_tms = pd.DataFrame(
                np.zeros((len(genes), len(celltypes) * len(ages_tms)), np.float32),
                index=genes,
                columns=columns_age,
                )
        ncells_exp_tp_tms = pd.Series(
                np.zeros(len(columns_age), np.int64), index=columns_age,
                )
        for ct in celltypes:
            print(ct)
            if ct not in celltypes_tms:
                # FIXME: set this as -1 or something?
                continue

            adata_ct = adata_tissue[adata_tissue.obs['cellType'] == ct]

            for age in ages_tms:
                idx_age = (adata_ct.obs['age'] == age).values.nonzero()[0]
                Xct_age = adata_ct.X[idx_age]
                label = '_'.join([ct, 'TMS', age])
                avg_exp_tp_tms[label] = np.asarray(Xct_age.mean(axis=0))[0]
                frac_exp_tp_tms[label] = np.asarray((Xct_age > 0).mean(axis=0))[0]
                ncells_exp_tp_tms[label] = len(idx_age)

    if switches['RNA-combine']:
        print('Combine Emily and TMS RNA data')
        def sorter(key):
            ct, dataset, tp = key.split('_')
            return (celltypes.index(ct), int(tp.strip('~m')), ['Emily', 'TMS'].index(dataset))

        # NOTE: using only Emily's genes as a "left" join
        avg_exp_tp_both = pd.merge(
                left=avg_exp_tp, right=avg_exp_tp_tms, how='left',
                left_index=True, right_index=True,
                )
        frac_exp_tp_both = pd.merge(
                left=frac_exp_tp, right=frac_exp_tp_tms, how='left',
                left_index=True, right_index=True,
                )
        ncells_exp_tp_both = pd.concat([ncells_exp_tp, ncells_exp_tp_tms])

        col_order = sorted(avg_exp_tp_both.columns.tolist(), key=sorter)
        avg_exp_tp_both = avg_exp_tp_both[col_order]
        frac_exp_tp_both = frac_exp_tp_both[col_order]
        ncells_exp_tp_both = ncells_exp_tp_both[col_order]

    if switches['RNA-addfeatures']:
        print('Add gene locations, will be useful to identify relationships with chromatin regions')
        genes = avg_exp.index.values

        import gzip
        with gzip.open('../../data/gene_annotations/mm10.ncbiRefSeq.gtf.gz', 'rt') as gtf:
            gene_annos = []
            for line in gtf:
                if '\ttranscript\t' not in line:
                    continue
                fields = line.split('\t')
                if fields[2] != 'transcript':
                    continue
                attrs = fields[-1].split(';')
                gene_name = None
                transcript_id = None
                for attr in attrs:
                    if 'gene_name' in attr:
                        gene_name = attr.split(' ')[-1][1:-1]
                    elif 'transcript_id' in attr:
                        transcript_id = attr.split(' ')[-1][1:-1]
                if (gene_name is None) or (transcript_id is None):
                    continue
                gene_annos.append({
                    'transcript_id': transcript_id,
                    'gene_name': gene_name,
                    'chromosome_name': fields[0],
                    'start_position': int(fields[3]),
                    'end_position': int(fields[4]),
                    'strand': 1 if fields[6] == '+' else -1,
                    'transcription_start_site': int(fields[3]) if fields[6] == '+' else int(fields[4]),
                    })
        gene_annos = pd.DataFrame(gene_annos)

        assert gene_annos['transcript_id'].value_counts()[0] == 1

        gene_annos.set_index('transcript_id', inplace=True)

        genes_missing = set(genes) - set(gene_annos['gene_name'].values)
        gene_annos_miss = pd.DataFrame([], index=genes_missing)
        gene_annos_miss['transcript_id'] = gene_annos_miss.index
        gene_annos_miss['start_position'] = -1
        gene_annos_miss['end_position'] = -1
        gene_annos_miss['strand'] = 0
        gene_annos_miss['chromosome_name'] = ''
        gene_annos_miss['transcription_start_site'] = -1
        gene_annos = pd.concat([gene_annos, gene_annos_miss])
        gene_annos = gene_annos.loc[genes]
        gene_annos['strand'] = gene_annos['strand'].astype('i2')

    if switches['RNA-store']:
        print('Store combined RNA data to file')
        with h5py.File(fn_out, 'a') as h5_data:
            ge = h5_data.create_group('gene_expression')
            ge.create_dataset('features', data=avg_exp.index.values.astype('S'))

            group = ge.create_group('feature_annotations')
            group.create_dataset('transcript_id', data=gene_annos.index.values.astype('S'))
            group.create_dataset('transcription_start_site', data=gene_annos['transcription_start_site'].values, dtype='i8')
            group.create_dataset('chromosome_name', data=gene_annos['chromosome_name'].astype('S'))
            group.create_dataset('start_position', data=gene_annos['start_position'].values, dtype='i8')
            group.create_dataset('end_position', data=gene_annos['end_position'].values, dtype='i8')
            group.create_dataset('strand', data=gene_annos['strand'].values, dtype='i2')

            group = ge.create_group('celltype')
            group.create_dataset('index', data=avg_exp.columns.values.astype('S'))
            group.create_dataset('average', data=avg_exp.T.values, dtype='f4')
            group.create_dataset('fraction', data=frac_exp.T.values, dtype='f4')
            group.create_dataset('cell_count', data=ncells_exp.values, dtype='i8')

            group = ge.create_group('celltype_dataset_timepoint')
            group.create_dataset('index', data=avg_exp_tp_both.columns.values.astype('S'))
            group.create_dataset('average', data=avg_exp_tp_both.T.values, dtype='f4')
            group.create_dataset('fraction', data=frac_exp_tp_both.T.values, dtype='f4')
            group.create_dataset('cell_count', data=ncells_exp_tp_both.values, dtype='i8')

            group = ge.create_group('celltype_dataset_timepoint_disease')
            group.create_dataset('index', data=avg_exp_di.columns.values.astype('S'))
            group.create_dataset('average', data=avg_exp_di.T.values, dtype='f4')
            group.create_dataset('fraction', data=frac_exp_di.T.values, dtype='f4')
            group.create_dataset('cell_count', data=ncells_exp_di.values, dtype='i8')
