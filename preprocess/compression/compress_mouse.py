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
import yaml
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import anndata
import scanpy as sc


tms_source_url = 'https://figshare.com/articles/dataset/Processed_files_to_use_with_scanpy_/8273102'
emilys_data_folder = pathlib.Path(__file__).parent.parent.parent / 'data' / 'Emilyslab'
tms_data_folder = pathlib.Path(__file__).parent.parent.parent / 'data' / 'tabula_muris_senis'
webapp_fdn = pathlib.Path('../../webapp')
output_fdn = webapp_fdn / 'static' / 'scData'
fn_out = output_fdn / 'mouse_compressed_heart_atlas.h5'


conversion_dict = {
    'TMS': {
        'endothelial cell of coronary artery': 'Coronary_EC',
        'fibroblast of cardiac tissue': 'Fibroblast',
        'endocardial cell': 'Endocardial',
        'smooth muscle cell': 'Smooth_muscle',
        'cardiac neuron': 'Neuronal',
        'mast cell': 'Myeloid',
        # See separate chat about vent VS atrial cardiomyocytes in TMS/Emily's
        'cardiomyocyte': 'Vent',
        'leukocyte': 'Leukocyte',
        'erythrocyte': 'Erythrocyte',
    },
    'Emily': {
        'Neuronal': 'Neuron',
        'Coronary_EC': 'Coronary',
        'Lymphatic_EC': 'Lymphatic',
        'Myeloid': 'Macrophage',
        'Bcell': 'B cell',
        'Tcell': 'T cell',
        'Smooth_muscle': 'Smooth muscle',
    }
}


def correct_celltypes(celltypes_raw, dataset):
    celltypes_new = celltypes_raw.copy()
    for ctraw, ct in conversion_dict[dataset].items():
        celltypes_new[celltypes_new == ctraw] = ct
    return celltypes_new



if __name__ == '__main__':

    print('Read celltype order from config file')
    with open(webapp_fdn / 'config.yml') as f:
        config = yaml.safe_load(f)
    celltypes = config['order']['celltypes']

    # Remove existing file if present
    if os.path.isfile(fn_out):
        os.remove(fn_out)

    switches = {
        'ATAC': True,
        'RNA': True,
        'RNA-TMS': True,
        'RNA-combine': True,
        'RNA-addfeatures': True,
        'RNA-store': True,
    }

    if switches['ATAC']:
        print('Emily\'s data (ATAC)')
        fn_atac = emilys_data_folder / 'Mouse_Heart_ATAC_peaks.h5ad'

        # raw access for speed/memory reasons
        with h5py.File(fn_atac) as fa:

            print('Read and check metadata')
            regions = fa['var']['_index'].asstr()[:]

            annos_a = fa['obs']['cellType'].asstr()[:]
            # FIX: ATAC-Seq and RNA-Seq both have a typo in lymphatics
            annos_a[annos_a == 'Lympathic_EC'] = "Lymphatic_EC"

            annos_a = correct_celltypes(annos_a, 'Emily')

            annos_a_count = pd.Series(annos_a).value_counts()

            celltypes_a = annos_a_count[~annos_a_count.index.str.startswith('unknown')].index

            assert set(celltypes) == set(celltypes_a)

            ages_a = fa['obs']['age'].asstr()[:]
            # Convert Y -> 10m, MA -> 19m
            ages_a[ages_a == 'Y'] = '10m'
            ages_a[ages_a == 'MA'] = '~19m'
            ages = ['10m', '~19m']

            assert(set(ages_a) == set(ages))

            # NOTE: "condition" is broken, reconstruct it from "group"
            condition = fa['obs']['group'].asstr()[:]
            condition[condition == 'MAE'] = 'E'
            condition[condition == 'YE'] = 'E'
            condition[condition == 'MAS'] = 'S'
            condition[condition == 'YS'] = 'S'

            columns_age = []
            for ct in celltypes:
                for age in ages:
                    columns_age.append('_'.join([ct, 'Emily', age]))

            print('Read ATAC matrix')
            Xh5 = fa['raw']['X']
            X = csr_matrix(
                (np.asarray(Xh5['data']), np.asarray(Xh5['indices']), np.asarray(Xh5['indptr'])),
                )

            print('Compress ATAC data for healthy cell types, time, and disease+time')
            avg_acc = pd.DataFrame(
                    np.zeros((len(regions), len(celltypes)), np.float32),
                    index=regions,
                    columns=celltypes)
            ncells_acc = pd.Series(np.zeros(len(celltypes), np.int64), index=celltypes)
            avg_acc_tp = pd.DataFrame(
                    np.zeros((len(regions), len(celltypes) * len(ages)), np.float32),
                    index=regions,
                    columns=columns_age,
                    )
            ncells_acc_tp = pd.Series(
                    np.zeros(len(columns_age), np.int64), index=columns_age,
                    )
            avg_acc_di = pd.DataFrame(
                    np.zeros((len(regions), len(celltypes) * len(ages)), np.float32),
                    index=regions,
                    columns=columns_age,
                    )
            ncells_acc_di = pd.Series(
                    np.zeros(len(columns_age), np.int64), index=columns_age,
                    )
            for ct in celltypes:
                print(ct)
                # Focus on baseline (sedentary) for now
                idx = ((annos_a == ct) & (condition == "S")).nonzero()[0]
                Xct = X[idx].astype(np.float32)

                # Binarize
                Xct.data[:] = 1

                # TODO: normalize
                adata_tmp = anndata.AnnData(X=Xct)
                sc.pp.normalize_total(
                    adata_tmp,
                    target_sum=None, # ?? Normalizing to the median cell coverage, does not make a lot of sense but also no harm
                    key_added='coverage',
                )
                Xct = adata_tmp.X

                # Avg accessibility == fraction accessible
                avg_acc[ct] = np.asarray(Xct.mean(axis=0))[0]

                # Number of cells
                ncells_acc[ct] = len(idx)

                # Split by time
                for age in ages:
                    idx_age = (ages_a[idx] == age).nonzero()[0]
                    Xct_age = Xct[idx_age]
                    label = '_'.join([ct, 'Emily', age])
                    avg_acc_tp[label] = np.asarray(Xct_age.mean(axis=0))[0]
                    ncells_acc_tp[label] = len(idx_age)

                # Disease/non-baseline (exercise)
                idx = ((annos_a == ct) & (condition == "E")).nonzero()[0]
                Xct = X[idx].astype(np.float32)
                for age in ages:
                    idx_age = (ages_a[idx] == age).nonzero()[0]
                    Xct_age = Xct[idx_age]
                    label = '_'.join([ct, 'Emily', age])
                    avg_acc_di[label] = np.asarray(Xct_age.mean(axis=0))[0]
                    ncells_acc_di[label] = len(idx_age)

        del X

        print('Store compressed atlas (ATAC)')
        with h5py.File(fn_out, 'a') as h5_data:
            ca = h5_data.create_group('chromatin_accessibility')
            ca.create_dataset('features', data=avg_acc.index.values.astype('S'))
            group = ca.create_group('celltype')
            group.create_dataset('index', data=avg_acc.columns.values.astype('S'))
            group.create_dataset('average', data=avg_acc.T.values, dtype='f4')
            group.create_dataset('cell_count', data=ncells_acc.values, dtype='i8')

            group = ca.create_group('celltype_dataset_timepoint')
            group.create_dataset('index', data=avg_acc_tp.columns.values.astype('S'))
            group.create_dataset('average', data=avg_acc_tp.T.values, dtype='f4')
            group.create_dataset('cell_count', data=ncells_acc_tp.values, dtype='i8')

            group = ca.create_group('celltype_dataset_timepoint_disease')
            group.create_dataset('index', data=avg_acc_di.columns.values.astype('S'))
            group.create_dataset('average', data=avg_acc_di.T.values, dtype='f4')
            group.create_dataset('cell_count', data=ncells_acc_di.values, dtype='i8')

    if switches['RNA']:
        print('Emily\'s data (RNA)')
        fn_rna = emilys_data_folder / 'Mouse_Heart_RNA.h5ad'

        # raw access for speed/memory reasons
        with h5py.File(fn_rna) as fr:

            print('Read and check metadata')
            genes = fr['var']['_index'].asstr()[:]

            annos_r = fr['obs']['cellType'].asstr()[:]
            # FIX: RNA-Seq uses "Pericytes" instead of "Pericyte"
            annos_r[annos_r == 'Pericytes'] = "Pericyte"
            # FIX: RNA-Seq has a typo in lymphatics
            annos_r[annos_r == 'Lympathic_EC'] = "Lymphatic_EC"

            annos_r = correct_celltypes(annos_r, 'Emily')

            annos_r_count = pd.Series(annos_r).value_counts()

            celltypes_r = annos_r_count[~annos_r_count.index.str.startswith('unknown')].index

            assert set(celltypes) == set(celltypes_r)

            ages_r = fr['obs']['age'].asstr()[:]
            # Convert Y -> 10m, MA -> 19m
            ages_r[ages_r == 'Y'] = '10m'
            ages_r[ages_r == 'MA'] = '~19m'
            ages = ['10m', '~19m']

            assert(set(ages_r) == set(ages))

            # NOTE: "condition" seems to be correct in the RNA data
            condition = fr['obs']['condition'].asstr()[:]

            columns_age = []
            for ct in celltypes:
                for age in ages:
                    columns_age.append('_'.join([ct, 'Emily', age]))

            print('Read RNA matrix')
            Xh5 = fr['raw']['X']
            X = csr_matrix(
                (np.asarray(Xh5['data']), np.asarray(Xh5['indices']), np.asarray(Xh5['indptr'])),
                )
            # FIXME: even though it says "raw", it's cptt -> natural logpp.
            # So we just undo the log and average (arithmetic after normalisation)
            X.data = np.exp(X.data) - 1

            print('Compress RNA data for healthy cell types, time, and disease+time')
            print('Do not store to file yet, we need to interleave it with TMS')
            avg_exp = pd.DataFrame(
                    np.zeros((len(genes), len(celltypes)), np.float32),
                    index=genes,
                    columns=celltypes)
            frac_exp = pd.DataFrame(
                    np.zeros((len(genes), len(celltypes)), np.float32),
                    index=genes,
                    columns=celltypes)
            ncells_exp = pd.Series(np.zeros(len(celltypes), np.int64), index=celltypes)
            avg_exp_tp = pd.DataFrame(
                    np.zeros((len(genes), len(celltypes) * len(ages)), np.float32),
                    index=genes,
                    columns=columns_age,
                    )
            frac_exp_tp = pd.DataFrame(
                    np.zeros((len(genes), len(celltypes) * len(ages)), np.float32),
                    index=genes,
                    columns=columns_age,
                    )
            ncells_exp_tp = pd.Series(
                    np.zeros(len(columns_age), np.int64), index=columns_age,
                    )
            avg_exp_di = pd.DataFrame(
                    np.zeros((len(genes), len(celltypes) * len(ages)), np.float32),
                    index=genes,
                    columns=columns_age,
                    )
            frac_exp_di = pd.DataFrame(
                    np.zeros((len(genes), len(celltypes) * len(ages)), np.float32),
                    index=genes,
                    columns=columns_age,
                    )
            ncells_exp_di = pd.Series(
                    np.zeros(len(columns_age), np.int64), index=columns_age,
                    )
            for ct in celltypes:
                print(ct)
                # Focus on baseline (sedentary) for now
                idx = ((annos_r == ct) & (condition == "S")).nonzero()[0]
                Xct = X[idx].astype(np.float32)

                # Avg expression
                avg_exp[ct].iloc[:Xct.shape[1]] = np.asarray(Xct.mean(axis=0))[0]

                # Fraction accessible
                frac_exp[ct].iloc[:Xct.shape[1]] = np.asarray((Xct > 0).mean(axis=0))[0]

                # Number of cells
                ncells_exp[ct] = len(idx)

                # Split by time
                for age in ages:
                    idx_age = (ages_r[idx] == age).nonzero()[0]
                    Xct_age = Xct[idx_age]
                    label = '_'.join([ct, 'Emily', age])
                    avg_exp_tp[label].iloc[:Xct_age.shape[1]] = np.asarray(Xct_age.mean(axis=0))[0]
                    frac_exp_tp[label].iloc[:Xct_age.shape[1]] = np.asarray((Xct_age > 0).mean(axis=0))[0]
                    ncells_exp_tp[label] = len(idx_age)

                # Disease/non-baseline (exercise)
                idx = ((annos_r == ct) & (condition == "E")).nonzero()[0]
                Xct = X[idx].astype(np.float32)
                for age in ages:
                    idx_age = (ages_r[idx] == age).nonzero()[0]
                    Xct_age = Xct[idx_age]
                    label = '_'.join([ct, 'Emily', age])
                    avg_exp_di[label].iloc[:Xct_age.shape[1]] = np.asarray(Xct_age.mean(axis=0))[0]
                    frac_exp_di[label].iloc[:Xct_age.shape[1]] = np.asarray((Xct_age > 0).mean(axis=0))[0]
                    ncells_exp_di[label] = len(idx_age)
            del X

    if switches['RNA-TMS']:
        print('Tabula Muris Senis')
        # Locate h5ad files for TMS
        drop_source = tms_data_folder / 'tabula-muris-senis-droplet-processed-official-annotations-Heart_and_Aorta.h5ad'

        # Ignore the FACS data for this repo: Emily used droplets anyway
        data_type = 'drop'
        data_filename = drop_source

        # Read the data
        adata_tissue = anndata.read(data_filename)

        # Restart from raw data and renormalize
        adata_tissue = adata_tissue.raw.to_adata()

        # cptt, like Emily's RNA data
        sc.pp.normalize_total(
            adata_tissue,
            target_sum=1e4,
            key_added='coverage',
        )

        #standardizing obs categories, cleaning up so that categories are just raw annotations
        adata_tissue.obs['cellType'] = correct_celltypes(
                np.asarray(adata_tissue.obs['cell_ontology_class']), 'TMS')
        celltypes_tms = adata_tissue.obs['cellType'].value_counts().index
        adata_tissue.obs['condition'] = 'S'

        # TODO: call northstar or treasuremap or both: for now, just take at face
        # value and rename some cell types to be somewhat consistent

        print('Add TMS data in celltype-timepoint group only')
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
