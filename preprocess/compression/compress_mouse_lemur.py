# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress mouse lemur atlas, heart.
'''
import os
import sys
import numpy as np
import pandas as pd
import h5py
import anndata


data_fdn = '../webapp/static/scData/'


def correct_annotations(adata, column):
    '''Overwrite annotations for macrophages and DCs'''
    annotations = pd.read_csv(
        '../data/tabula_microcebus/reannotations.tsv',
        sep='\t',
        index_col=0,
    ).squeeze(axis=1)
    annotations = pd.Series(annotations, dtype='category')

    cats_old = adata.obs[column].cat.categories
    cats_new = list(set(annotations.cat.categories) - set(cats_old))
    adata.obs[column] = adata.obs[column].cat.add_categories(cats_new)
    annotations = annotations.cat.set_categories(adata.obs[column].cat.categories)

    adata.obs.loc[annotations.index, column] = annotations
    adata.obs[column] = adata.obs[column].cat.remove_unused_categories()


if __name__ == '__main__':

    # Load data
    print('Load single cell data')
    fn_atlas = '../data/tabula_microcebus/Heart_FIRM_hvg.h5ad'
    adata = anndata.read_h5ad(fn_atlas)

    # NOTE: the lemur data is in some weird normalization between 0 and 8.79,
    # and adata.raw is None. Hmm, there must be a logp1 there, let's try to
    # undo that transformation by hand. After np.expm1 the sum of each cell is
    # 10000, so there you go, multiply by 100 and done.
    adata.X = 100 * np.expm1(adata.X)

    print('Compress')
    # Find cell type column
    column = 'cell_ontology_class_v1'
    # NOTE: there is a free annotation column that is informative, basically
    # a higher resolution clustering but with some questionable things like
    # CD4+ CD8+ T cells. Let's go trad for now

    print('Reannotate aerocytes, DCs, etc.')
    correct_annotations(adata, column)

    # Exclude unassigned/doublets
    unwanted_types = [
        'lymphocyte',
        'unassigned',
        'epithelial cell of uterus',
        'mature NK T cell',
        'endothelial cell',
        'granulocyte monocyte progenitor cell',
    ]
    adata = adata[~adata.obs[column].isin(unwanted_types)]

    # Set cell type order
    #celltypes = adata.obs[column].unique()
    celltypes = [
        'Mesothelial',
        'Adventitial FB',
        'Alveolar FB',
        'ASM',
        'VSM',
        'Pericyte',
        'Lymphatic',
        'Arterial II',
        'Venous',
        'gCap',
        'Aerocyte',
        'Schwann cell',
        'Plasmablast',
        'B cell',
        'NK cell',
        'T cell',
        'IL cell',
        'pDC',
        'DC I',
        'DC II',
        'DC III',
        'Alveolar mac',
        'Interstitial mac',
        'Monocyte',
        'Neutrophil',
        'Basophil',
        'Eosinophil',
        'Platelet',
        'Alveolar type I',
        'Alveolar type II',
        'Ciliated',
        'Club',
        'Brush',
        'Basal',
    ]

    # Merge some annotations together
    annotation_merges = {
        'Mesothelial': ['mesothelial cell'],
        'Adventitial FB': ['fibroblast'],
        'Alveolar FB': ['fibroblast of lung'],
        # NOTE: they express HHIP, no PDGFRA...
        'ASM': ['myofibroblast cell'],
        'VSM': ['vascular associated smooth muscle cell'],
        'Pericyte': ['pericyte cell'],
        'Schwann cell': [
            'myelinating Schwann cell',
            'non-myelinating Schwann cell',
        ],
        'Venous': ['vein endothelial cell'],
        'gCap': ['capillary endothelial cell'],
        'Aerocyte': ['Aerocyte'],
        'Arterial II': ['endothelial cell of artery'],
        'Lymphatic': ['endothelial cell of lymphatic vessel'],
        'Plasmablast': ['plasma cell'],
        'NK cell': ['natural killer cell'],
        'IL cell': ['innate lymphoid cell'],
        'B cell': ['B cell'],
        'T cell': [
            'T cell',
            'CD4-positive, alpha-beta T cell',
            'CD8-positive, alpha-beta T cell',
        ],
        'red blood cell lineage': [
            'erythroid lineage cell',
            'hematopoietic precursor cell',
            'erythroid progenitor cell',
            'megakaryocyte progenitor cell',
        ],
        'Alveolar type I': ['type I pneumocyte'],
        'Alveolar type II': ['type II pneumocyte'],
        'Ciliated': ['lung ciliated cell'],
        'Club': ['club cell'],
        'Brush': ['brush cell of bronchus'],
        'Basal': ['basal cell of epithelium of bronchus'],
        'DC I': ['DC I'],
        'DC II': ['DC II'],
        'DC III': ['DC III'],
        'pDC': ['plasmacytoid dendritic cell'],
        'Alveolar mac': ['alveolar macrophage'],
        'Interstitial mac': ['macrophage'],
        'Monocyte': ['monocyte'],
        'Neutrophil': ['neutrophil'],
        'Basophil': ['basophil'],
        'Eosinophil': ['eosinophil'],
        'Platelet': ['platelet'],
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
        # Avg gene expression (compute cptt after averaging)
        avg_exp[ct] = np.asarray(adata[adata.obs[column+'_compressed'] == ct].X.mean(axis=0))[0]
        avg_exp[ct] *= 1e4 / avg_exp[ct].sum()
        # Proportion expressing
        frac_exp[ct] = np.asarray((adata[adata.obs[column+'_compressed'] == ct].X > 0).mean(axis=0))[0]
        # Number of cells
        ncells[ct] = (adata.obs[column+'_compressed'] == ct).sum()

    # Store to file
    print('Store compressed atlas')
    fn_out = data_fdn + 'mouselemur_compressed_heart_atlas.h5'
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
    output_fn = data_fdn + 'human_gene_friends.h5'
    with h5py.File(output_fn, 'w') as f:
        for gene, friends_i in friends.items():
            lmax = max(len(x) for x in friends_i)
            friends_i = friends_i.astype('S'+str(lmax))
            dset = f.create_dataset(gene, data=friends_i)

