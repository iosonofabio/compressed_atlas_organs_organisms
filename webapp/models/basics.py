import h5py
import numpy as np
import pandas as pd

from models.assets import fn_atlasd


def read_tissues(
        species='mouse',
        feature_type='gene_expression',
        ):
    '''Read the list of tissues from the atlas file'''
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        tissues = np.array(h5_data[feature_type]['tissues'])
    return tissues


def read_feature_order(
        feature_type='gene_expression',
        species='mouse',
        ):
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        features = np.array(h5_data[feature_type]["features"].asstr())
    feature_order = pd.Series(data=np.arange(len(features)), index=features)
    return feature_order


def read_cell_types(
        tissue,
        feature_type='gene_expression',
        species='mouse',
        ):
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        celltypes = np.array(
            h5_data[feature_type]["by_tissue"][tissue]['celltype']['index'].asstr(),
            )

    return np.asarray(celltypes)


def read_gene_order(species='mouse'):
    return read_feature_order(
            feature_type='gene_expression',
            species=species,
            )


def read_regions_order(species='mouse'):
    return read_feature_order(
            feature_type='chromatin_accessibility',
            species=species,
            )
