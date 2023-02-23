# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/04/22
content:    Some quirks of the celltypes in the data.
'''
import numpy as np
from config import configuration as config

from models.basics import (
    read_tissues,
    read_cell_types,
)


def validate_correct_celltypestr(
        celltypestr,
        species='mouse',
    ):
    '''Validate cell type names and correct misspellings if possible'''

    # TODO: check punctuation more accurately
    celltypes = celltypestr.strip(' ').replace('.', ',').replace(';', ',').split(',')
    
    celltypesv = {}
    tissues = read_tissues(species=species)
    for tissue in tissues:
        celltypes_tissue = read_cell_types(
                tissue, species=species,
            )

        for celltype in celltypes:
            if celltype in celltypes_tissue:
                celltypesv[celltype] = celltype
                continue

            if (celltype.endswith('s')) and (celltype in celltypes_tissue):
                celltypesv[celltype] = celltype[:-1]
                continue

        if len(celltypesv) == len(celltypes):
            break

    else:
        # Some cell types did not validate

        # TODO: implement spelling correction
        return None

    celltypestr = ','.join(celltypesv[x] for x in celltypes)

    return celltypestr


def validate_correct_celltypedatasettimepoint(suffix):
    '''Validate celltype + dataset + timepoint'''
    #TODO
    return suffix

