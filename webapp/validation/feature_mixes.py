# vim: fdm=indent
'''
author:     Fabio Zanini
date:       24/12/22
content:    Validate mixes of different types of features, e.g. gene expression
            and chromatin accessibility.
'''
import numpy as np
import pandas as pd

from validation.genes import validate_correct_gene
from validation.regions import validate_correct_multiregion



def validate_correct_feature_mix(featurestr, species='mouse'):
    '''Validate and correct a mix features.

    For now, we just use heuristics to guess which features come from which
    type of assay. This is robust enough for text recognition but probably not
    for voice control.

    featurestr: A string containing the features.
    species: Which species do the feature belong to.
    '''
    features_raw = featurestr.replace(' ', '').split(',')

    features_df = pd.DataFrame([], index=features_raw)
    features_df['raw'] = features_df.index
    features_df['index'] = np.arange(len(features_raw))
    ftypes = []
    for feature_raw in features_raw:
        if '-' in feature_raw:
            ft = 'chromatin_accessibility'
        elif ':' in feature_raw:
            ft = 'chromatin_accessibility'
        else:
            ft = 'gene_expression'
        ftypes.append(ft)
    features_df['type_raw'] = ftypes
    features_df['type'] = ''
    features_df['corrected'] = ''

    for feature_raw, type_raw in features_df['type_raw'].items():
        if type_raw == 'gene_expression':
            featurec = validate_correct_gene(feature_raw, species=species)
            if featurec is None:
                continue
            features_df.at[feature_raw, 'corrected'] = featurec
            features_df.at[feature_raw, 'type'] = type_raw
        elif type_raw == 'chromatin_accessibility':
            # Multiregion can return more than one region (e.g. an interval),
            # however they are all of the same feature type
            featurec = validate_correct_multiregion(feature_raw, species=species)
            if featurec is None:
                continue
            features_df.at[feature_raw, 'corrected'] = featurec
            features_df.at[feature_raw, 'type'] = type_raw
        else:
            raise ValueError('Heuristic for unimplemented feature type')

    features_df = features_df.loc[features_df['type'] != '']

    result = {}
    for feature_type, group in features_df.groupby('type'):
        # This is needed because of multiregion
        features_cand = (','.join(group['corrected'].values.astype(str))).split(',')
        # Avoid duplicates
        features = []
        for fea in features_cand:
            if fea not in features:
                features.append(fea)
        result[feature_type] = ','.join(features)

    return result

