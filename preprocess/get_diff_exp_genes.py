# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/04/22
content:    Precompute list of top correlated genes for "friends" button.
'''
import os
import sys
import h5py
import numpy as np
import pandas as pd
from scipy import stats


data_fdn = '../webapp/static/scData/'
fn_atlas = data_fdn + 'condensed_lung_atlas_ordered.h5'


def get_counts(condition='normal'):
    category = 'celltype_dataset_timepoint'
    if condition == 'hyperoxia':
        category += '_hyperoxia'

    with h5py.File(fn_atlas) as f:
        dic = f[category]['gene_expression_average']
        genes = np.array(dic['axis0'].asstr())
        celltypes = np.array(dic['axis1'].asstr())
        counts = np.array(dic['block0_values']).astype(np.float32)

    # Fix this gene that has "inf" counts: CT010467.1
    counts[:, 5370] = 0
    counts = 1e6 * (counts.T / counts.sum(axis=1)).T

    # FIXME: focus on highly expressed genes
    idx = counts.max(axis=0) > 30
    counts = counts[:, idx]
    genes = genes[idx]
    counts = pd.DataFrame(counts.T, index=genes, columns=celltypes)

    return counts


def align_count_tables(counts, counts_ho):
    # FIXME: ho has fewer genes because of a bug
    overlap = list(set(counts_ho.index) & set(counts.index))
    counts = counts.loc[overlap]
    counts_ho = counts_ho.loc[overlap]

    # Pad columns with zeros
    cols = list(set(counts_ho.columns) | set(counts.columns))
    for col in cols:
        if col not in counts_ho.columns:
            counts_ho[col] = 0
        if col not in counts.columns:
            counts[col] = 0
    counts = counts[cols]
    counts_ho = counts_ho[cols]

    return counts, counts_ho


if __name__ == '__main__':


    # First comparison: hyperoxia/normal
    counts = get_counts()
    counts_ho = get_counts(condition='hyperoxia')

    counts, counts_ho = align_count_tables(counts, counts_ho)

    print('Compute log2FC')
    log2fc = np.log2(counts_ho + 0.5) - np.log2(counts + 0.5)

    print('Compute diff exp')
    nct = log2fc.shape[1]
    diff_exp = {}
    for i, col in enumerate(log2fc.columns):
        print(f'Column {i+1}/{nct}', end='\r')
        top_up = log2fc[col].nlargest(20).index.tolist()
        top_down = log2fc[col].nsmallest(20).index.tolist()
        diff_exp[col] = {
            'up': top_up,
            'down': top_down,
            'both': top_up + top_down,
        }
    print()

    print('Store to file')
    output_fn = data_fdn + 'degs.h5'
    with h5py.File(output_fn, 'w') as f:
        group = f.create_group('hyperoxia')
        for i, ct in enumerate(diff_exp):
            groupi = group.create_group(ct)
            for key in ['up', 'down', 'both']:
                genes_i = np.asarray(diff_exp[ct][key])
                lmax = max(len(x) for x in genes_i)
                genes_i = genes_i.astype('S'+str(lmax))
                groupi.create_dataset(key, data=genes_i)
