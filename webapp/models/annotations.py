import h5py
import numpy as np
import pandas as pd

from models.assets import fn_atlasd


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
        regions = np.array(h5_data['chromatin_accessibility']["features"].asstr())
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
