# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/04/22
content:    Some quirks of the celltypes in the data.
'''
import numpy as np
from config import configuration as config


# A few globals with cell type manipulation data structures
celltypes_validated = list(config['order']['celltypes'])

celltype_tuples = list(celltypes_validated)
if ('conversions' in config) and ('celltypes' in config['conversions']):
    for celltype_sloppy, celltype in config['conversions']['celltypes'].items():
        idx = celltype_tuples.index(celltype)
        celltype_tuples[idx] = (celltype_sloppy, celltype)

celltype_dict = {}
for ct in celltype_tuples:
    if isinstance(ct, str):
        celltype_dict[ct] = ct
    else:
        celltype_dict[ct[0]] = ct[1]
celltype_dict_inv = {val: key for key, val in celltype_dict.items()}


def adjust_celltypes(celltypes_raw, species='mouse'):
    # TODO: reorder even when it's stratified by dataset/timepoint
    for ct in celltypes_raw:
        if '_Emily' in ct:
            return celltypes_raw, np.arange(len(celltypes_raw))

    celltypes_raw = list(celltypes_raw)
    ct_adj = []
    idx = []
    for ct in celltype_tuples:
        if isinstance(ct, str):
            ct1, ct2 = ct, ct
        else:
            ct1, ct2 = ct
        if ct2 in celltypes_raw:
            idx.append(celltypes_raw.index(ct2))
            ct_adj.append(ct2)
        elif ct1 in celltypes_raw:
            idx.append(celltypes_raw.index(ct1))
            ct_adj.append(ct2)
    return np.array(ct_adj), np.array(idx)


def rename_celltypes(celltypes_raw, inverse=False):
    '''Rename celltypes according to the standard table above, no order change'''
    if not inverse:
        return [celltype_dict[x] for x in celltypes_raw]
    else:
        return [celltype_dict_inv[x] for x in celltypes_raw]


def validate_correct_celltypestr(celltypestr):
    '''Validate cell type names and correct misspellings if possible'''
    # TODO: check punctuation more accurately
    celltypes = celltypestr.strip(' ').replace('.', ',').replace(';', ',').split(',')

    # Capitalization is not great, lowercase check
    celltypes = [x.lower() for x in celltypes]
    celltypesd = {}
    for ct1, ct2 in celltype_dict.items():
        celltypesd[ct1.lower()] = ct2
        celltypesd[ct2.lower()] = ct2

    # Validate
    celltypesv = []
    for celltype in celltypes:
        if celltype in celltypesd:
            celltypesv.append(celltypesd[celltype])
            continue

        # Cut plural if it is found
        if celltype.endswith('s'):
            celltype = celltype[:-1]
            if celltype in celltypesd:
                celltypesv.append(celltypesd[celltype])
                continue

        # TODO: implement more spelling correction
        return None

    celltypestr = ','.join(celltypesv)
    return celltypestr


def validate_correct_celltypedatasettimepoint(suffix):
    '''Validate celltype + dataset + timepoint'''
    #TODO
    return suffix

