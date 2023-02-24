# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Interpret text from Google TTS into a redirect URL.
'''
from urllib.parse import urlencode
from flask import url_for

from models import (
        get_features_nearby,
        get_marker_features,
        get_differential_url,
        )


command_dict = {
    'expression_by_celltype': {
        'url_func': lambda sfx: url_for(
            'measurement_by_celltype')+'?featurestring='+sfx,
    },
    'accessibility_by_celltype': {
        'url_func': lambda sfx: url_for(
            'measurement_by_celltype')+'?featurestring='+sfx,
    },
    'expression_overtime_1feature': {
        'url_func': lambda sfx: url_for(
            'measurement_overtime_1feature')+'?featurestring='+sfx,
    },
    'accessibility_overtime_1feature': {
        'url_func': lambda sfx: url_for(
            'measurement_overtime_1feature')+'?featurestring='+sfx,
    },
    'gene_friends': {
        'url_func': 'TODO',  # TODO
    },
    'nearby_regions': {
        'url_func': lambda sfx, distance=0: url_for('measurement_by_celltype') + '?' + \
                'featurestring='+get_features_nearby(
                    sfx,
                    target_types=['chromatin_accessibility'],
                    dmax=distance,
                    include_query=True,
                    ),
    },
    'nearby_genes': {
        'url_func': lambda sfx, distance=0: url_for('measurement_by_celltype') + '?' + \
                'featurestring='+get_features_nearby(
                    sfx,
                    target_types=['gene_expression'],
                    dmax=distance,
                    include_query=True,
                    ),
    },
    'nearby_features': {
        'url_func': lambda sfx, distance=0: url_for('measurement_by_celltype') + '?' + \
                'featurestring='+get_features_nearby(
                    sfx,
                    target_types=['chromatin_accessibility', 'gene_expression'],
                    dmax=distance,
                    include_query=True,
                    ),
    },
    'marker_genes': {
        'url_func': lambda sfx, number=10: url_for('measurement_by_celltype') + '?' + \
                'featurestring='+get_marker_features(
                    sfx,
                    feature_type="gene_expression",
                    ntop=number,
                    ),
    },
    'marker_regions': {
        'url_func': lambda sfx, number=10: url_for('measurement_by_celltype') + '?' + \
                'featurestring='+get_marker_features(
                    sfx,
                    feature_type="chromatin_accessibility",
                    ntop=number,
                    ),
    },
    'differentially_expressed_genes': {
        'url_func': lambda sfx, number=20: url_for(
            'heatmap_differential_measurement')+"?"+get_differential_url(
                sfx, feature_type='gene_expression', kind='both',
                n_features=number)
    },
    'differentially_accessible_regions': {
        'url_func': lambda sfx, number=20: url_for(
            'heatmap_differential_measurement')+"?"+get_differential_url(
                sfx, feature_type='chromatin_accessibility', kind='both',
                n_features=number)
    },
    'differentially_measured_features': {
        'url_func': lambda sfx, number=20: url_for(
            'heatmap_differential_measurement')+"?"+get_differential_url(
                sfx, feature_type='all', kind='both', n_features=number)
    },
    'upregulated_genes': {
        'url_func': lambda sfx, number=10: url_for(
            'heatmap_differential_measurement')+"?"+get_differential_url(
                sfx, feature_type='gene_expression', kind='up',
                n_features=number)
    },
    'upregulated_regions': {
        'url_func': lambda sfx, number=10: url_for(
            'heatmap_differential_measurement')+"?"+get_differential_url(
                sfx, feature_type='chromatin_accessibility', kind='up',
                n_features=number)
    },
    'upregulated_features': {
        'url_func': lambda sfx, number=10: url_for(
            'heatmap_differential_measurement')+"?"+get_differential_url(
                sfx, feature_type='all', kind='up',
                n_features=number)
    },
    'downregulated_genes': {
        'url_func': lambda sfx, number=10: url_for(
            'heatmap_differential_measurement')+"?"+get_differential_url(
                sfx, feature_type='gene_expression', kind='down',
                n_features=number)
    },
    'downregulated_regions': {
        'url_func': lambda sfx, number=10: url_for(
            'heatmap_differential_measurement')+"?"+get_differential_url(
                sfx, feature_type='chromatin_accessibility', kind='down',
                n_features=number)
    },
    'downregulated_features': {
        'url_func': lambda sfx, number=10: url_for(
            'heatmap_differential_measurement')+"?"+get_differential_url(
                sfx, feature_type='all', kind='down',
                n_features=number)
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
    '''Add organismal species to url if needed'''
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

    # Some categories benefit from optional additional parameters, e.g. the
    # number of DEGs, the distance from a gene
    kwargs = {}
    if ('marker' in category) or ('differential' in category):
        if 'number' in text_dict['prefix_parameters']:
            kwargs['number'] = text_dict['prefix_parameters']['number']

    if 'nearby' in category:
        kwargs['distance'] = text_dict['prefix_parameters'].get('distance', 50000)

    url = command_dict[category]['url_func'](suffix_corrected, **kwargs)

    if 'species=' not in url:
        url = add_species_to_url(url, text_dict['species'])

    return {
        'outcome': 'success',
        'url': url,
    }
