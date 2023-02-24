# vim: fdm=indent
'''
author:     Fabio Zanini
date:       09/05/22
content:    Interpret text about differential expression.
'''
import numpy as np

from config import configuration as config
from .celltypes import (
        validate_correct_celltypestr,
        )


def sanitize_differential_text_after_vs(text_after_vs):
    '''Sanitize text after "vs" to fix typical typos'''
    condition = config['condition'].lower()
    baseline = config.get('baseline', 'normal').lower()

    # Autocorrect conditions, they are usually long words
    for condword in [condition, baseline]:
        if len(text_after_vs) == len(condword):
            char1 = np.array(list(text_after_vs))
            char2 = np.array(list(condword))
            if (char1 != char2).sum() <= 2:
                return condword

    return text_after_vs


def sanitize_differential_suffix(suffix):
    '''Sanitize suffix for differential expression'''
    condition = config['condition'].lower()
    baseline = config.get('baseline', 'normal').lower()

    suffix = (suffix.replace(' VS ', ' vs ')
                    .replace(config['condition'], condition)
                    .replace(config.get('baseline', 'normal'), baseline))

    return suffix


def none_in(substrings, string):
    '''Check whether none of the substrings are found in the string'''
    for substring in substrings:
        if substring in string:
            return False
    return True


def _parse_suffix_no_versus(suffix):
    '''Parse suffix when no VS clause is found'''
    condition = config['condition'].lower()
    baseline = config.get('baseline', 'normal').lower()

    # Default to mouse, TODO: enable other species here
    species = list(config['paths']['compressed_atlas'].keys())[0]
    datasets_all = config['order']['dataset'][species]

    # If found, cut "in XXX (condition)" from the suffix
    suffix = suffix.replace('in '+condition, '')

    # Find time point
    deli_tp = ' at '
    if deli_tp not in suffix:
        # Assume P7 ACZ for the short syntax
        tmp = config['defaults']['condition']['dataset_timepoint'][0]
        dataset, tp = tmp.split('_')
        prefix_nodataset = suffix
    else:
        idx_deli_tp = suffix.find(deli_tp)
        tp = suffix[idx_deli_tp+len(deli_tp):].split(' ')[0]
        # Find data set
        for dataset in datasets_all:
            if ' in '+dataset in suffix:
                deli_ds = ' in '+dataset
                idx_deli_ds = suffix.find(deli_ds)
                prefix_nodataset = suffix[:idx_deli_ds] + suffix[idx_deli_ds+len(deli_ds):]
                break
        else:
            # Guess first dataset from timepoints
            # TODO: this could improve
            dataset = datasets_all[0]
            prefix_nodataset = suffix

    # Find out cell type
    deli_ct = ' in '
    idx_deli_ct = prefix_nodataset.find(deli_ct)
    if idx_deli_ct == -1:
        raise ValueError('Sentence structure for DEG not recognized')
    celltype_raw = prefix_nodataset[idx_deli_ct+len(deli_ct):]
    celltype = validate_correct_celltypestr(celltype_raw)

    conditions = [
        {'celltype': celltype,
         'dataset': dataset,
         'timepoint': tp,
         'condition': condition,
         },
        {'celltype': celltype,
         'dataset': dataset,
         'timepoint': tp,
         'conditino': baseline,
         },
    ]
    return {
        'kind': 'condition/baseline',
        'conditions': conditions,
    }


def _parse_suffix_timepoints(suffix, idx_deli, text_after_vs):
    '''Parse suffix for differential measurement between 2 time points'''
    species = list(config['paths']['compressed_atlas'].keys())[0]
    timepoints_all = config['order']['timepoint'][species]
    datasets_all = config['order']['dataset'][species]

    tp2 = text_after_vs
    # The word before vs should be the first time point
    tp1 = suffix[:idx_deli].split()[-1]
    if tp1 not in timepoints_all:
        raise ValueError(
            'Sentence structure for differential measurement not recognized')
    prefix_notp = suffix[:suffix.find(' at ')]

    # Find out dataset
    for dataset in datasets_all:
        if ' in '+dataset in suffix[:idx_deli]:
            deli_ds = ' in '+dataset
            idx_deli_ds = prefix_notp.find(deli_ds)
            prefix_nodataset = prefix_notp[:idx_deli_ds] + \
                    prefix_notp[idx_deli_ds+len(deli_ds):idx_deli]
            break
    else:
        # Guess first dataset from timepoints
        dataset = datasets_all[0]
        prefix_nodataset = prefix_notp

    # Find out cell type
    deli_ct = ' in '
    idx_deli_ct = prefix_nodataset.find(deli_ct)
    if idx_deli_ct == -1:
        raise ValueError(
            'Sentence structure for differential measurement not recognized')
    celltype_raw = prefix_nodataset[idx_deli_ct+len(deli_ct):]
    celltype = validate_correct_celltypestr(celltype_raw)

    conditions = [
        {'celltype': celltype,
         'dataset': dataset,
         'timepoint': tp1,
         },
        {'celltype': celltype,
         'dataset': dataset,
         'timepoint': tp2,
         },
    ]
    return {
        'kind': 'timepoints',
        'conditions': conditions,
    }


def _parse_suffix_datasets(suffix, idx_deli, text_after_vs):
    raise NotImplementedError("DEG between data sets not implemented")


def _parse_suffix_condition_baseline(suffix, idx_deli, text_after_vs):
    '''Parse suffix for differential measurement condition vs baseline'''
    condition = config['condition'].lower()
    baseline = config.get('baseline', 'normal').lower()

    species = list(config['paths']['compressed_atlas'].keys())[0]
    timepoints_all = config['order']['timepoint'][species]
    datasets_all = config['order']['dataset'][species]

    # Find time point
    deli_tp = ' at '
    if deli_tp not in suffix:
        # Assume P7 ACZ for the short syntax
        tmp = config['defaults']['condition']['dataset_timepoint'][0]
        dataset, tp = tmp.split('_')
        prefix_nodataset = suffix
    else:
        idx_deli_tp = suffix.find(deli_tp)
        tp = suffix[idx_deli_tp+len(deli_tp):].split(' ')[0]
        # Find data set
        for dataset in datasets_all:
            if ' in '+dataset in suffix:
                deli_ds = ' in '+dataset
                idx_deli_ds = suffix.find(deli_ds)
                prefix_nodataset = suffix[:idx_deli_ds] + \
                    suffix[idx_deli_ds+len(deli_ds):]
                break
        else:
            # Guess first dataset from timepoints
            dataset = datasets_all[0]
            prefix_nodataset = suffix

    # Find out cell type
    deli_ct = ' in '
    idx_deli_ct = prefix_nodataset.find(deli_ct)
    if idx_deli_ct == -1:
        raise ValueError('Sentence structure for DEG not recognized')
    celltype_raw = prefix_nodataset[idx_deli_ct+len(deli_ct):]
    idx_endct1 = celltype_raw.find(' in ')
    idx_endct2 = celltype_raw.find(' at ')
    if (idx_endct1 == -1) and (idx_endct2 == -1):
        pass
    elif idx_endct1 == -1:
        celltype_raw = celltype_raw[:idx_endct2]
    elif idx_endct2 == -1:
        celltype_raw = celltype_raw[:idx_endct1]
    else:
        celltype_raw = celltype_raw[:min(idx_endct1, idx_endct2)]

    celltype = validate_correct_celltypestr(celltype_raw)

    conditions = [
        {'celltype': celltype,
         'dataset': dataset,
         'timepoint': tp,
         'condition': condition,
         },
        {'celltype': celltype,
         'dataset': dataset,
         'timepoint': tp,
         'condition': baseline,
         },
    ]
    if condition in text_after_vs:
        conditions = conditions[::-1]

    return {
        'kind': 'condition/baseline',
        'conditions': conditions,
    }


def _parse_suffix_celltypes(suffix, idx_deli, text_after_vs):
    '''Parse suffix for differential measurement between 2 cell types'''
    deli2 = ' in '
    idx_deli2 = suffix[:idx_deli].rfind(deli2)
    celltypestr = suffix[idx_deli2+len(deli2):idx_deli]+','+text_after_vs
    celltypestr_vali = validate_correct_celltypestr(celltypestr)
    celltypes = celltypestr_vali.split(',')
    if len(celltypes) != 2:
            raise ValueError(
                    'Sentence structure for differential measurement not recognized')
    ct1, ct2 = celltypes
    # Find time point
    deli_tp = ' at '
    if deli_tp not in suffix[:idx_deli]:
        # No time point, just averages
        conditions = [{'celltype': ct1}, {'celltype': ct2}]
        return conditions
    idx_deli_tp = suffix[:idx_deli].find(deli_tp)
    tp = suffix[idx_deli_tp+len(deli_tp):idx_deli].split(' ')[0]
    # Find data set
    for dataset in datasets_all:
        if ' in '+dataset in suffix[:idx_deli]:
            break
    else:
        # Guess first dataset
        dataset = datasets_all[0]

    conditions = [
        {'celltype': ct1,
         'dataset': dataset,
         'timepoint': tp,
         },
        {'celltype': ct2,
         'dataset': dataset,
         'timepoint': tp,
         },
    ]
    return {
        'kind': 'celltypes',
        'conditions': conditions,
    }


def get_differential_conditions(suffix):
    '''Interpret the differential measurement conditions from a suffix'''
    condition = config['condition'].lower()
    baseline = config.get('baseline', 'normal').lower()

    suffix = sanitize_differential_suffix(suffix)

    if none_in((' vs ', ' versus ', condition), suffix):
        raise ValueError('Conditions not recognized')

    # Default to mouse, TODO: enable other species here
    species = list(config['paths']['compressed_atlas'].keys())[0]
    timepoints_all = config['order']['timepoint'][species]

    # Default to condition vs baseline for a cell type and (default) timepoint
    # and dataset
    if none_in((' vs ', ' versus '), suffix):
        return _parse_suffix_no_versus(suffix)

    # Figure out dataset, timepoint, and cell type
    # The sentence structure is such that whatever is after "vs" determines
    # the actual comparison. For now we assume that all stays constant except
    # for one thing. So let's start there
    deli = ' vs '
    if deli not in suffix:
        deli = ' versus '
    idx_deli = suffix.find(deli)
    text_after_vs = suffix[idx_deli+len(deli):]

    # Autocorrect typical sentence endings
    text_after_vs = sanitize_differential_text_after_vs(text_after_vs)

    # Comparison between time points
    if text_after_vs in timepoints_all:
        return _parse_suffix_timepoints(suffix, idx_deli, text_after_vs)

    # Comparison between datasets
    if text_after_vs in config['order']['dataset']:
        return _parse_suffix_datasets(suffix, idx_deli, text_after_vs)

    # Condition/baseline
    if text_after_vs in (condition, baseline):
        return _parse_suffix_condition_baseline(
                suffix, idx_deli, text_after_vs)

    # Two cell types
    return _parse_suffix_celltypes(suffix, idx_deli, text_after_vs)
