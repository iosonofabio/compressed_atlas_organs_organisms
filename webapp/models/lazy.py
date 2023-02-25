import h5py
import numpy as np
import pandas as pd

from config import configuration as config
from models.assets import fn_atlasd

feature_types = config['feature_types']

# A bunch of lazy data that can be read from file only when needed
data_lazy = {
    'speciess': None,
    'tissuesd': None,
    'celltypesd': None,
    'feature_orderd': None,
    'feature_annotationd': None,
    'gene_matrixd': None,
}


def read_speciess():
    '''Read list of species'''
    return list(fn_atlasd.keys())


def load_speciess():
    '''Preload species'''
    data_lazy['speciess'] = read_speciess()


def read_tissues(
        feature_type='gene_expression',
        species='mouse',
        ):
    '''Read list of tissues from the atlas file for a species and modality'''
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        tissues = np.array(h5_data[feature_type]['tissues'].asstr())
    return tissues


def load_tissues():
    '''Preload tissues for each species and modality'''
    tissuesd = {}
    for ftype in feature_types:
        tissuesd[ftype] = {}
        for species in fn_atlasd:
            tissuesd[ftype][species] = read_tissues(
                    feature_type=ftype,
                    species=species,
                    )
    data_lazy['tissuesd'] = tissuesd


def read_cell_types(
        tissue,
        feature_type='gene_expression',
        species='mouse',
        ):
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        if tissue is not None:
            celltypes = np.array(
                h5_data[feature_type]["by_tissue"][tissue]['celltype']['index'].asstr(),
                )
        else:
            supertypes = h5_data[feature_type]['celltypes']['supertypes'].asstr()
            celltypes = []
            for supertype in supertypes:
                tmp = h5_data[feature_type]['celltypes'][supertype].asstr()
                celltypes.extend(list(tmp))

    return np.asarray(celltypes)


def load_celltypes():
    '''Preload cell types for each tissue and species'''

    celltypesd = {}
    for ftype in feature_types:
        celltypesd[ftype] = {}
        for species in fn_atlasd:
            tissues = get_tissues(ftype, species)
            celltypesd[ftype][species] = {}
            for tissue in list(tissues) + [None]:
                celltypesd[ftype][species][tissue] = read_cell_types(
                    tissue,
                    feature_type=ftype,
                    species=species,
                    )
    data_lazy['celltypesd'] = celltypesd


def read_feature_order(
        feature_type='gene_expression',
        species='mouse',
        ):
    '''Read feature order from file'''
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        features = np.array(h5_data[feature_type]["features"].asstr())
    feature_order = pd.Series(data=np.arange(len(features)), index=features)
    return feature_order


def load_feature_orderd():
    feature_orderd = {}
    for ftype in feature_types:
        feature_orderd[ftype] = {}
        for species in fn_atlasd:
            feature_orderd[ftype][species] = read_feature_order(
                    species=species,
                    feature_type=ftype,
                    )
    data_lazy['feature_orderd'] = feature_orderd


def read_gene_annotations(species='mouse'):
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        genes = np.array(h5_data['gene_expression']["features"].asstr())
        group = h5_data['gene_expression']['feature_annotations']
        start = group['start_position'][:]
        end = group['end_position'][:]
        chromosome = np.array(group['chromosome_name'].asstr())
        strand = group['strand'][:]
        tss = group['transcription_start_site'][:]
    annotations = pd.DataFrame({
        'start': start, 'end': end, 'chrom': chromosome, 'strand': strand,
        'tss': tss,
    }, index=genes)
    return annotations


def read_region_annotations(species='mouse'):
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        regions = np.array(
                h5_data['chromatin_accessibility']["features"].asstr())
    df = pd.DataFrame([], index=regions)
    tmp = df.index.str.split('-', expand=True)
    df['chrom'] = tmp.get_level_values(0)
    df['start'] = tmp.get_level_values(1).astype(int)
    df['end'] = tmp.get_level_values(2).astype(int)
    # Use integer division to speed up computations a little bit
    df['mid'] = (df['start'] + df['end']) // 2
    return df


def read_feature_annotations(feature_type, species='mouse'):
    if feature_type == 'gene_expression':
        return read_gene_annotations(species=species)
    elif feature_type == 'chromatin_accessibility':
        return read_region_annotations(species=species)
    else:
        raise ValueError(f'Feature type not found: {feature_type}')


def load_feature_annotationd():
    feature_annotationd = {}
    for ft in feature_types:
        feature_annotationd[ft] = {}
        for species in fn_atlasd:
            feature_annotationd[ft][species] = read_feature_annotations(
                    feature_type=ft,
                    species=species,
                    )
    data_lazy['feature_annotationd'] = feature_annotationd


def load_gene_matrixd():
    '''Load the character matrix table for each gene and species

    This is used to autocorrect gene names.
    '''
    gene_orderd = get_feature_orderd('gene_expression')
    gene_matrixd = {}
    for species in gene_orderd:
        mtx = (pd.Series(gene_orderd[species].index)
                 .str
                 .split('', expand=True)
                 .iloc[:, 1:-1]
                 .fillna('')
                 .values
                 .astype('U1'))
        gene_matrixd[species] = mtx
    data_lazy['gene_matrixd'] = gene_matrixd


#####################################
# Getters
#####################################
def get_speciess():
    '''Return cached species'''
    if data_lazy['speciess'] is None:
        load_speciess()
    return data_lazy['speciess']


def get_tissues(feature_type, species):
    '''Return cached tissues'''
    if data_lazy['tissuesd'] is None:
        load_tissues()
    return data_lazy['tissuesd'][feature_type][species]


def get_celltypes(feature_type, species, tissue=None):
    '''Return cached celltypes'''
    if data_lazy['celltypesd'] is None:
        load_celltypes()
    return data_lazy['celltypesd'][feature_type][species][tissue]


def get_feature_orderd(feature_type, species=None):
    '''Return cached feature order dictionary'''
    if data_lazy['feature_orderd'] is None:
        load_feature_orderd()
    if species is None:
        return data_lazy['feature_orderd'][feature_type]
    return data_lazy['feature_orderd'][feature_type][species]


def get_feature_annotationd(feature_type, species=None):
    if data_lazy['feature_annotationd'] is None:
        load_feature_annotationd()
    if species is None:
        return data_lazy['feature_annotationd'][feature_type]
    return data_lazy['feature_annotationd'][feature_type][species]


def get_gene_matrixd(species):
    if data_lazy['gene_matrixd'] is None:
        load_gene_matrixd()
    return data_lazy['gene_matrixd'][species]


def get_feature_annotation(feature, species=None):
    '''Get annotation for a single feature'''
    for feature_type in feature_types:
        feature_annotationd = get_feature_annotationd(
            feature_type, species=species,
        )
        if feature in feature_annotationd.index:
            return {
                'feature_type': feature_type,
                'annotation': feature_annotationd.loc[feature],
                }

    return None


def get_gene_annotationd(species):
    return get_feature_annotation('gene_expression', species)


def get_region_annotationd(species):
    return get_feature_annotation('chromatin_accessibility', species)
