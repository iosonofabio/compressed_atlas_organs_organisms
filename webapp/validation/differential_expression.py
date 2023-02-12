# vim: fdm=indent
'''
author:     Fabio Zanini
date:       09/05/22
content:    Interpret text about differential expression.
'''
from .celltypes import (
        validate_correct_celltypestr,
        rename_celltypes,
        )
from config import configuration as config


def get_deg_conditions(suffix):
    '''Interpret the deg conditions from a suffix'''
    if ((' vs ' not in suffix) and
        (' versus ' not in suffix) and
        ('hyperoxia' not in suffix)):
        raise ValueError('Conditions not recognized')

    # Default to mouse, TODO: enable other species here
    timepoints_all = config['order']['timepoints']['mouse']

    # Default to hyperoxia vs normoxia for a specific time point,
    # cell type, and dataset
    if (' vs ' not in suffix) and (' versus ' not in suffix):
        # Cut "in hyperoxia" from the suffix
        suffix = suffix.replace('in hyperoxia', '')

        # Find time point
        deli_tp = ' at '
        if deli_tp not in suffix:
            # Assume P7 ACZ for the short syntax
            tp = 'P7'
            dataset = 'ACZ'
            prefix_nodataset = suffix
        else:
            idx_deli_tp = suffix.find(deli_tp)
            tp = suffix[idx_deli_tp+len(deli_tp):].split(' ')[0]
            # Find data set
            for dataset in ['ACZ', 'Hurskainen2021', 'TMS']:
                if ' in '+dataset in suffix:
                    deli_ds = ' in '+dataset
                    idx_deli_ds = suffix.find(deli_ds)
                    prefix_nodataset = suffix[:idx_deli_ds] + \
                            suffix[idx_deli_ds+len(deli_ds):]
                    break
            else:
                # Guess dataset from timepoints
                if tp.endswith('m'):
                    dataset = 'TMS'
                elif tp in ['E18.5', 'P1', 'P21']:
                    dataset = 'ACZ'
                else:
                    dataset = 'Hurskainen2021'
                prefix_nodataset = suffix

        # Find out cell type
        deli_ct = ' in '
        idx_deli_ct = prefix_nodataset.find(deli_ct)
        if idx_deli_ct == -1:
            raise ValueError('Sentence structure for DEG not recognized')
        celltype_raw = prefix_nodataset[idx_deli_ct+len(deli_ct):]
        celltype = rename_celltypes(
                [validate_correct_celltypestr(celltype_raw)],
                inverse=True)[0]

        conditions = [
            {'celltype': celltype, 'dataset': dataset,
             'timepoint': tp, 'disease': 'hyperoxia'},
            {'celltype': celltype, 'dataset': dataset,
             'timepoint': tp, 'disease': 'normal'},
        ]
        return {'kind': 'hyperoxia/normal', 'conditions': conditions}

    # Figure out dataset, timepoint, and cell type
    # The sentence structure is such that whatever is after "vs" determines
    # the actual comparison. For now we assume that all stays constant except
    # for one thing. So let's start there
    deli = ' vs '
    if deli not in suffix:
        deli = ' versus '
    idx_deli = suffix.find(deli)
    text_after_vs = suffix[idx_deli+len(deli):]
    if text_after_vs in timepoints_all:
        tp2 = text_after_vs
        # The word before vs should be the first time point
        tp1 = suffix[:idx_deli].split()[-1]
        if tp1 not in timepoints_all:
            raise ValueError('Sentence structure for DEG not recognized')
        prefix_notp = suffix[:suffix.find(' at ')]

        # Find out dataset
        for dataset in ['ACZ', 'Hurskainen2021', 'TMS']:
            if ' in '+dataset in suffix[:idx_deli]:
                deli_ds = ' in '+dataset
                idx_deli_ds = prefix_notp.find(deli_ds)
                prefix_nodataset = prefix_notp[:idx_deli_ds] + \
                        prefix_notp[idx_deli_ds+len(deli_ds):idx_deli]
                break
        else:
            # Guess dataset from timepoints
            prefix_nodataset = prefix_notp
            if tp1.endswith('m') or tp2.endswith('m'):
                dataset = 'TMS'
            elif (tp1 in ['E18.5', 'P1', 'P21']) or (tp2 in ['E18.5', 'P1', 'P21']):
                dataset = 'ACZ'
            else:
                dataset = 'Hurskainen2021'

        # Find out cell type
        deli_ct = ' in '
        idx_deli_ct = prefix_nodataset.find(deli_ct)
        if idx_deli_ct == -1:
            raise ValueError('Sentence structure for DEG not recognized')
        celltype_raw = prefix_nodataset[idx_deli_ct+len(deli_ct):]
        print(celltype_raw)
        celltype = rename_celltypes(
            [validate_correct_celltypestr(celltype_raw)], inverse=True)[0]

        conditions = [
            {'celltype': celltype, 'dataset': dataset, 'timepoint': tp1},
            {'celltype': celltype, 'dataset': dataset, 'timepoint': tp2},
        ]
        return {'kind': 'timepoints', 'conditions': conditions}

    # Comparison between datasets
    if text_after_vs in ['ACZ', 'TMS', 'Hurskainen2021']:
        raise NotImplementedError("DEG between data sets not implemented")

    # Hyperoxia/normoxia
    if text_after_vs in ('hyperoxia', 'normal'):
        # Find time point
        deli_tp = ' at '
        if deli_tp not in suffix:
            # Assume P7 ACZ for the short syntax
            tp = 'P7'
            dataset = 'ACZ'
            prefix_nodataset = suffix
        else:
            idx_deli_tp = suffix.find(deli_tp)
            tp = suffix[idx_deli_tp+len(deli_tp):].split(' ')[0]
            # Find data set
            for dataset in ['ACZ', 'Hurskainen2021', 'TMS']:
                if ' in '+dataset in suffix:
                    deli_ds = ' in '+dataset
                    idx_deli_ds = suffix.find(deli_ds)
                    prefix_nodataset = suffix[:idx_deli_ds] + \
                            suffix[idx_deli_ds+len(deli_ds):]
                    break
            else:
                # Guess dataset from timepoints
                if tp.endswith('m'):
                    dataset = 'TMS'
                elif tp in ['E18.5', 'P1', 'P21']:
                    dataset = 'ACZ'
                else:
                    dataset = 'Hurskainen2021'
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
        elif (idx_endct1 == -1):
            celltype_raw = celltype_raw[:idx_endct2]
        elif (idx_endct2 == -1):
            celltype_raw = celltype_raw[:idx_endct1]
        else:
            celltype_raw = celltype_raw[:min(idx_endct1, idx_endct2)]
        celltype = rename_celltypes(
            [validate_correct_celltypestr(celltype_raw)], inverse=True)[0]

        conditions = [
            {'celltype': celltype, 'dataset': dataset,
             'timepoint': tp, 'disease': 'hyperoxia'},
            {'celltype': celltype, 'dataset': dataset,
             'timepoint': tp, 'disease': 'normal'},
        ]
        if 'hyperoxia' in text_after_vs:
            conditions = conditions[::-1]
        return {'kind': 'hyperoxia/normal', 'conditions': conditions}

    # Two cell types
    deli2 = ' in '
    idx_deli2 = suffix[:idx_deli].rfind(deli2)
    celltypestr = suffix[idx_deli2+len(deli2):idx_deli]+','+text_after_vs
    print(celltypestr)
    celltypestr_vali = validate_correct_celltypestr(celltypestr)
    celltypes = celltypestr_vali.split(',')
    if len(celltypes) != 2:
            raise ValueError('Sentence structure for DEG not recognized')
    celltypes = rename_celltypes(celltypes, inverse=True)
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
    for dataset in ['ACZ', 'Hurskainen2021', 'TMS']:
        if ' in '+dataset in suffix[:idx_deli]:
            break
    else:
        # Guess dataset from timepoints
        if tp.endswith('m'):
            dataset = 'TMS'
        elif tp in ['E18.5', 'P1', 'P21']:
            dataset = 'ACZ'
        else:
            dataset = 'Hurskainen2021'
    conditions = [
        {'celltype': ct1, 'dataset': dataset, 'timepoint': tp},
        {'celltype': ct2, 'dataset': dataset, 'timepoint': tp},
    ]
    return {'kind': 'celltypes', 'conditions': conditions}
