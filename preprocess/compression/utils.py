'''
Utility functions for the compression
'''
import os
import gzip
import numpy as np
import pandas as pd
import h5py

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
        # TMC
        if tissue.endswith('_hvg'):
            tissue = tissue[:-len('_FIRM_hvg')].replace('_', ' ').title()

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
        ('human', 'immune cell'): {
            'T': ['CD3D', 'CD3G', 'CD3E', 'TRAC', 'IL7R'],
            'B': ['MS4A1', 'CD19', 'CD79A'],
            'NK': ['PFN1', 'TMSB4XP8'],
            'macrophage': ['MRC1', 'MARCO', 'CD163', 'C1QA', 'C1QB', 'CST3'],
            'dendritic': ['FCER1A', 'IL1R2', 'CD86', 'HLA-DPB1', 'HLA-DRB1'],
            'neutrophil': ['S100A8', 'S100A7'],
        },
        ('human', 'endothelial'): {
            'arterial': ['GJA5', 'BMX', 'SEMA3G', 'VIM'],
            'venous': ['VWF', 'MMRN2', 'CLEC14A', 'ACKR1'],
            'lymphatic': ['LYVE1', 'PROX1', 'THY1', 'MMRN1', 'TFF3', 'TFPI'],
            'capillary': ['SLC9A3R2', 'PLPP1', 'PECAM1', 'IGKC', 'CALD1', 'CRHBP', 'KDR'],
            'epithelial': ['COBLL1', 'EPCAM', 'CD24'],
            '': [
                'JUN', 'JUND', 'SQSTM1', 'SELENOH', 'FOS', 'ACP1', 'EPB41L2',
                'MALAT1', 'CAP1', 'FABP5P7', 'XIST', 'TGFBR2', 'SPARCL1',
                'FCN3', 'F8'],
            'acinar': ['PRSS2', 'ENPP2', 'GALNT15', 'APOD', 'CLPS'],
        },
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
        ('mouselemur', 'lymphocyte'): {
            'B': ['MS4A1', 'CD79A', 'CD79B', 'CD19'],
            'T': ['TRAC', 'CD3E', 'CD3D', 'CD3G'],
            'NK': ['GZMA', 'NCAM1', 'FCER1G', 'GZMK', 'KLRB1'],
            'macrophage': ['C1QA', 'CD68', 'MARCO', 'CST3'],
            'monocyte': ['PSAP', 'CD14'],
            'neutrophil': ['S100A8', 'S100A9', 'STFA1', 'STFA2'],
            'erythrocyte': ['BETA-S', 'ALAS2', 'HBB-B2', 'TMEM14C'],
            '': ['SNRPF'],
        },
    }

    bad_prefixes = {
        'mouse': ['Rpl', 'Rps', 'Linc', 'Mt'],
        'human': [
            'RPL', 'RPS', 'LINC', 'MT', 'EPAS1', 'DYNLL1',
            'EIF3G', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E',
            'GZMA', 'GNLY', 'CD74', 'KRT4', 'TYROBP'],
        'mouselemur': [
            'RPL', 'RPS', 'LINC', 'MT', 'EPAS1', 'DYNLL1',
            'EIF3G', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E',
            'GZMA', 'GNLY', 'CD74', 'KRT4', 'TYROBP',
            'UBA52', 'LOC1', 'MYBL2', 'MAL', 'ATP5A1', 'ARHGAP15'],
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
            found_bad_prefix = False
            for bad_pfx in bad_prefixes[species]:
                if gene.startswith(bad_pfx):
                    found_bad_prefix = True
                    break
            if found_bad_prefix:
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


def fix_annotations(adata, column, species, tissue, rename_dict, coarse_cell_types):
    '''Correct cell types in each tissue according to known dict'''
    blacklist = {
        ('mouselemur', 'Bone Marrow'): ['lymphocyte', 'type II pneumocyte'],
        ('mouselemur', 'Heart'): ['type II pneumocyte'],
        ('mouselemur', 'Kidney'): ['stromal cell', 'urothelial cell'],
        ('mouselemur', 'Lung'): ['epithelial cell of uterus'],
        ('mouselemur', 'Pancreas'): ['stromal cell', 'pancreatic endocrine cell'],
        ('mouselemur', 'Tongue'): ['stromal cell', 'pancreatic endocrine cell'],
    }

    celltypes_new = np.asarray(adata.obs[column]).copy()

    # Exclude blacklisted
    if (species, tissue) in blacklist:
        for ctraw in blacklist[(species, tissue)]:
            celltypes_new[celltypes_new == ctraw] = ''

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
                adata_coarse_type, species, celltype)

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


def collect_gene_annotations(anno_fn, genes):
    '''Collect gene annotations from GTF file'''
    with gzip.open(anno_fn, 'rt') as gtf:
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

    # FIXME: choose the largest transcript or something. For this particular
    # repo it's not that important
    gene_annos = (gene_annos.drop_duplicates('gene_name')
                            .set_index('gene_name', drop=False))

    genes_missing = list(set(genes) - set(gene_annos['gene_name'].values))
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

    return gene_annos


def store_compressed_atlas(
        fn_out,
        compressed_atlas,
        tissues,
        gene_annos,
        celltype_order,
        ):
    '''Store compressed atlas into h5 file'''
    genes = gene_annos.index.tolist()

    with h5py.File(fn_out, 'a') as h5_data:
        ge = h5_data.create_group('gene_expression')
        ge.create_dataset('features', data=np.array(genes).astype('S'))

        group = ge.create_group('feature_annotations')
        group.create_dataset(
                'gene_name', data=gene_annos.index.values.astype('S'))
        group.create_dataset(
                'transcription_start_site',
                data=gene_annos['transcription_start_site'].values, dtype='i8')
        group.create_dataset(
                'chromosome_name',
                data=gene_annos['chromosome_name'].astype('S'))
        group.create_dataset(
                'start_position',
                data=gene_annos['start_position'].values, dtype='i8')
        group.create_dataset(
                'end_position',
                data=gene_annos['end_position'].values, dtype='i8')
        group.create_dataset(
                'strand', data=gene_annos['strand'].values, dtype='i2')

        ge.create_dataset('tissues', data=np.array(tissues).astype('S'))

        supergroup = ge.create_group('by_tissue')
        for tissue in tissues:
            tgroup = supergroup.create_group(tissue)
            for label in ['celltype', 'celltype_dataset_timepoint']:
                avg_ge = compressed_atlas[tissue][label]['avg']
                frac_ge = compressed_atlas[tissue][label]['frac']
                ncells_ge = compressed_atlas[tissue][label]['ncells']

                group = tgroup.create_group(label)
                group.create_dataset(
                        'index', data=avg_ge.columns.values.astype('S'))
                group.create_dataset(
                        'average', data=avg_ge.T.values, dtype='f4')
                group.create_dataset(
                        'fraction', data=frac_ge.T.values, dtype='f4')
                group.create_dataset(
                        'cell_count', data=ncells_ge.values, dtype='i8')

        ct_group = ge.create_group('celltypes')
        supertypes = np.array([x[0] for x in celltype_order])
        ct_group.create_dataset(
                'supertypes',
                data=supertypes.astype('S'),
                )
        for supertype, subtypes in celltype_order:
            ct_group.create_dataset(
                supertype,
                data=np.array(subtypes).astype('S'),
            )
