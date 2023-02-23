'''
Utility functions for the compression
'''
import os
import numpy as np
import pandas as pd

import scanpy as sc


def get_tissue_data_dict(atlas_folder, rename_dict):
    '''Get a dictionary with tissue order and files'''
    result = []

    fns = os.listdir(atlas_folder)
    fns = [x for x in fns if '.h5ad' in x]

    for filename in fns:
        # TMS
        tissue = filename.split('-')[-1].split('.')[0]
        # TS
        if tissue.startswith('TS_'):
            tissue = tissue[3:]

        tissue = rename_dict['tissues'].get(tissue, tissue)
        result.append({
            'tissue': tissue,
            'filename': atlas_folder / filename,
        })

    result = pd.DataFrame(result).set_index('tissue')['filename']

    # Order tissues alphabetically
    result = result.sort_index()

    return result


def subannotate(adata, species, annotation, verbose=True):
    '''This function subannotates a coarse annotation from an atlasi.

    This is ad-hoc, but that's ok for now. Examples are 'lymphocyte', which is
    a useless annotation unless you know what kind of lymphocytes these are, or
    if it's a mixed bag.
    '''

    markers = {
        ('mouse', 'endothelial cell'): {
            'arterial': ['Gja5', 'Bmx'],
            'venous': ['Slc6a2', 'Vwf'],
            'lymphatic': ['Ccl21a', 'Prox1', 'Thy1'],
            'capillary': ['Rn45s', 'Slc6a6', 'Comt'],
            'smooth muscle': [
                'Thy1', 'Mustn1', 'Gng11', 'Mgp', 'Acta2', 'Aspn', 'Myl9'],
            'pericyte': ['Pdgfrb', 'Cox4i2', 'Higd1b'],
            'dendritic': ['Cd34', 'Cd300lg', 'Ly6c1', 'Ramp3'],
            'beta': ['Iapp', 'Ins1', 'Ins2', 'Srgn', 'Syngr2', 'Tsc22d1',
                     'Igfbp3'],
            'alpha': ['Chga', 'Gcg'],
            'acinar': ['Prss2', 'Try5', 'Sycn', 'Ctrb1', 'Clps', 'Ndrg1', 'Fabp4'],
            'stellate': ['Plac9'],
            'PP': ['Ppy'],
        },
        ('mouse', 'lymphocyte'): {
            'B': ['Ms4a1', 'Cd79a', 'Cd79b', 'Cd19'],
            'T': ['Trac', 'Cd3e', 'Cd3d', 'Cd3g'],
            'NK': ['Gzma', 'Ncam1'],
            'macrophage': ['C1qa', 'Cd68', 'Marco', 'Cst3'],
            'monocyte': ['Psap', 'Cd14'],
            'neutrophil': ['S100a8', 'S100a9', 'Stfa1', 'Stfa2'],
            'erythrocyte': ['Beta-s', 'Alas2', 'Hbb-b2', 'Tmem14c'],
            '': ['Snrpf'],
        },
        ('mouse', 'leukocyte'): {
            'B': ['Ms4a1', 'Cd79a', 'Cd79b', 'Cd19'],
            'T': ['Trac', 'Cd3e', 'Cd3d', 'Cd3g'],
            'NK': ['Gzma', 'Ncam1'],
            'macrophage': ['C1qa', 'Cd68', 'Marco', 'Cst3'],
            'monocyte': ['Psap', 'Cd14'],
            'neutrophil': ['S100a8', 'S100a9', 'Stfa1', 'Stfa2'],
        },
    }

    markersi = markers.get((species, annotation), None)
    if markersi is None:
        raise ValueError(
            f'Cannot subannotate without markers for {species}, {annotation}')

    adata = adata.copy()
    sc.pp.log1p(adata)

    genes, celltypes = [], []
    for celltype, markers_ct in markersi.items():
        celltypes.append(celltype)
        for gene in markers_ct:
            if gene in adata.var_names:
                genes.append(gene)
            elif verbose:
                print('Missing gene:', gene)

    adatam = adata[:, genes].copy()

    # No need for PCA because the number of genes is small

    # Get neighbors
    sc.pp.neighbors(adatam)

    # Get communities
    sc.tl.leiden(adatam)

    adata.obs['subleiden'] = adatam.obs['leiden']
    sc.tl.rank_genes_groups(adata, 'subleiden')
    top_marker = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(2)

    subannos = {}
    for cluster, genestop in top_marker.items():
        found = False
        for gene in genestop:
            if found:
                break
            if gene.startswith('Rps') or gene.startswith('Rpl'):
                subannos[cluster] = ''
                continue
            for celltype, markers_ct in markersi.items():
                if gene in markers_ct:
                    subannos[cluster] = celltype
                    found = True
                    break
            else:
                import ipdb; ipdb.set_trace()
                raise ValueError('Marker not found:', gene)
        if not found:
            subannos[cluster] = ''

    new_annotations = adata.obs['subleiden'].map(subannos)

    return new_annotations


def correct_celltypes(adata, column, species):
    '''Correct cell types in each tissue according to known dict'''
    celltypes_new = np.asarray(adata.obs[column]).copy()

    # Rename according to standard dict
    for ctraw, celltype in rename_dict['cell_types'].items():
        celltypes_new[celltypes_new == ctraw] = celltype

    ct_found = np.unique(celltypes_new)

    # Look for coarse annotations
    ctnew_list = set(celltypes_new)
    for celltype in ctnew_list:
        if celltype in coarse_cell_types:
            idx = celltypes_new == celltype
            adata_coarse_type = adata[idx]
            subannotations = subannotate(
                adata_coarse_type, 'mouse', celltype)

            # Ignore reclustering into already existing types, we have enough
            for subanno in subannotations:
                if subanno in ct_found:
                    subannotations[subannotations == subanno] = ''

            celltypes_new[idx] = subannotations

    return celltypes_new


def get_celltype_order(celltypes_unordered, celltype_order):
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
