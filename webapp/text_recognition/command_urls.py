# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Interpret text from Google TTS into a redirect URL.
'''
import numpy as np
import pandas as pd
from flask import url_for
from urllib.parse import urlencode

from models import get_marker_features, get_de_url



command_dict = {
    'expression_by_celltype': {
        'url_func': lambda sfx: url_for('measurement_by_celltype')+'?featurestring='+sfx,
    },
    'accessibility_by_celltype': {
        'url_func': lambda sfx: url_for('measurement_by_celltype')+'?featurestring='+sfx,
    },
    'expression_overtime_1feature': {
        'url_func': lambda sfx: url_for('expression_overtime_1feature')+'?featurestring='+sfx,
    },
    'gene_friends': {
        'url_func': 'TODO',  # TODO
    },
    'marker_genes': {
        'url_func': lambda sfx: url_for('measurement_by_celltype') + '?' + \
                'genestring='+get_marker_features(sfx),
    },
    'differentially_expressed_genes': {
        'url_func': lambda sfx: url_for('heatmap_differential_genes')+"?"+get_de_url(sfx, kind='both')
    },
    'upregulated_genes': {
        'url_func': lambda sfx: url_for('heatmap_differential_genes')+"?"+get_de_url(sfx, kind='up')
    },
    'downregulated_genes': {
        'url_func': lambda sfx: url_for('heatmap_differential_genes')+"?"+get_de_url(sfx, kind='down')
    },
    'list_cell_types': {
        'url_func': lambda sfx: url_for(
            'list_celltypes_timepoint',
            timepoint=sfx),
    },
    'celltype_abundance': {
        'url_func': lambda sfx: url_for(
            'plot_celltype_abundance',
            timepoint=sfx),
    },
    'compare_species': {
        'url_func': lambda dic: url_for(
            'heatmap_species_comparison') + '?' + urlencode(dic)
    },
}


def add_species_to_url(url, species):
    # Mouse is still the default
    if species == 'mouse':
        return url

    if '?' not in url:
        url += '?species='+species
    else:
        url += '&species='+species
    return url


def get_command_response(text_dict):
    '''Format response to each given command'''
    if text_dict is None:
        # Default answer raises an alert in the frontend
        return {'outcome': 'fail'}

    suffix_corrected = text_dict['suffix_corrected']
    question = text_dict['question']
    category = text_dict['category']

    if suffix_corrected is None:
        return {
            'outcome': 'question',
            'question': question,
            question: text_dict['suffix'],
            'url_prefix': command_dict[category]['url_func'](''),
        }

    url = command_dict[category]['url_func'](suffix_corrected)
    if 'species=' not in url:
        url = add_species_to_url(url, text_dict['species'])

    return {
        'outcome': 'success',
        'url': url,
    }

