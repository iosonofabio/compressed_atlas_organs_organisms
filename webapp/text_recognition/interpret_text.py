# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Interpret text from Google TTS into a redirect URL.
'''
from validation.genes import validate_correct_genestr
from validation.regions import (
    validate_correct_single_region,
    validate_correct_multiregion,
    )
from validation.celltypes import (
    validate_correct_celltypestr,
    validate_correct_celltypedatasettimepoint,
    )
from validation.timepoints import validate_correct_timepoint
from validation.species import (
    validate_correct_species,
    validate_correct_species_genestr,
)


phrase_dict = {
    'expression_by_celltype': {
        'prefixes': [
            'what is the expression of',
            'what\'s the expression of',
            'what is the gene expression of',
            'what\'s the gene expression of',
            'expression of',
            'gene expression of',
            'show expression of',
            'show gene expression of',
            'show the expression of',
            'show the gene expression of',
            'exp of',
            'e of',
            '!e',
        ],
        'suffix_type': 'genestring',
    },
    'accessibility_by_celltype': {
        'prefixes': [
            'what is the accessibility of',
            'what\'s the accessibility of',
            'what is the chromatin accessibility of',
            'accessibility of',
            'chromatin accessibility of',
            'show the chromatin accessibility of',
            'show chromatin accessibility of',
            'show accessibility of',
            'show the accessibility of',
            'acc of',
            'a of',
            '!a',
        ],
        'suffix_type': 'regionmultistring',
    },
    'expression_overtime_1feature': {
        'prefixes': [
            'what is the developmental expression of',
            'what\'s the developmental expression of',
            'what is the developmental progression of',
            'what\'s the developmental progression of',
            'what is the gene progression of',
            'what\'s the gene progression of',
            'progression of',
            'developmental expression of',
            'developmental progression of',
            'expression across development of',
            'show developmental progression of',
            'show developmental expression of',
            'show the developmental expression of',
            'show the developmental gene expression of',
            'show the progression of',
            'show progression of',
            'prog of',
            'p of',
            '!p',
        ],
        'suffix_type': 'genestring',
    },
    'gene_friends': {
        'prefixes': [
            'what are gene friends of',
            'gene friends of',
            'correlates of',
            'show the gene friends of',
            'show gene friends of',
            'show correlates of',
            'friends of',
            'f of',
            '!f',
        ],
        'suffix_type': 'genestring',
    },
    'marker_genes': {
        'prefixes': [
            'what are the markers of',
            'what are markers of',
            'what are marker genes of',
            'markers for',
            'markers of',
            'marker genes of',
            'marker genes for',
            'show markers of',
            'show markers for',
            'show marker genes of',
            'show marker genes for',
            'show the markers of',
            'show the markers for',
            'show the marker genes of',
            'show the marker genes for',
            '!m',
        ],
        'suffix_type': 'celltypestring',
    },
    'differentially_expressed_genes': {
        'prefixes': [
            'differentially expressed genes',
            'degs',
            '!d',
        ],
        'suffix_type': 'celltype_dataset_timepoint_string',
    },
    'upregulated_genes': {
        'prefixes': [
            'genes upregulated in',
            'upregulated genes in',
            'upregulated in',
            'up in',
            '!u',
        ],
        'suffix_type': 'celltype_dataset_timepoint_string',
    },
    'downregulated_genes': {
        'prefixes': [
            'genes downregulated in',
            'downregulated genes in',
            'downregulated in',
            'down in',
            '!d',
        ],
        'suffix_type': 'celltype_dataset_timepoint_string',
    },
    'list_cell_types': {
        'prefixes': [
            'list cell types',
            'list the cell types',
            'list cell type ',
            'list the cell type '
            'show cell types',
            'what are the cell types',
            'what cell types are there',
            'what cell types are present',
            'cell types',
            '!ct',
        ],
        'suffix_type': 'timepoint',
    },
    'celltype_abundance': {
        'prefixes': [
            'show cell type abundance at',
            'cell type abundance at',
            'abundance of cell types at',
        ],
        'suffix_type': 'timepoint',
    },
    'compare_species': {
        'prefixes': [
            'compare expression of',
            'compare gene expression of',
            'compare the expression of',
            'compare the gene expression of',
        ],
        'suffix_type': 'species_genestring',
    },
}
phrase_dict_inv = {}
for key, val in phrase_dict.items():
    for prefix in val['prefixes']:
        phrase_dict_inv[prefix] = key


def infer_command_from_text(text_raw):
    '''Figure out category of command'''
    # Cut punctuation at the end of the command
    text = text_raw.rstrip('?!.')

    # Prefixes are checked case-insensitive
    text_lower = text.strip().lower()

    for prefix, category in phrase_dict_inv.items():
        if not text_lower.startswith(prefix):
            continue

        # Figure out type of suffix
        suffix_type = phrase_dict[category]['suffix_type']

        # Cut prefix... there are two common cases for how to treat the suffix
        cats_keep_whitespace = (
            'celltype_dataset_timepoint_string',
            'timepoint',
            'celltypestring',
            'species_genestring',
            )
        if suffix_type in cats_keep_whitespace:
            suffix = text[len(prefix):]
        else:
            suffix = text_lower[len(prefix):]

        # Extract species from suffix
        suffix, species = excise_species_from_suffix(suffix)

        if suffix_type not in cats_keep_whitespace:
            # Remove whitespace from suffix
            suffix = suffix.replace(' ', '')

        return {
            'prefix': prefix,
            'suffix': suffix,
            'category': category,
            'species': species,
            }
    return None


def excise_phrase(text, phrase):
    '''Excise a phrase from a text'''
    if phrase not in text:
        return text
    if text.endswith(phrase):
        # Assume space before
        return text[:-len(phrase)-1]
    # Assume space after
    idx = text.find(phrase+' ')
    return text[:idx]+text[idx+len(phrase)+1:]


def excise_species_from_suffix(suffix):
    '''Excise species from suffix if found'''
    phrases = [
        # NOTE: order matters
        ('human', ['in humans', 'in human']),
        ('lemur', ['in mouse lemur', 'in lemur', 'in monkey']),
        ('mouse', ['in mouse', 'in mice']),
    ]
    for species, phrases_species in phrases:
        for phrase in phrases_species:
            if phrase in suffix:
                suffix = excise_phrase(suffix, phrase)
                return suffix, species

    # Default species is still mouse
    return suffix, 'mouse'


def interpret_text(text):
    '''Interpret natural language text as command'''
    inferred_dict = infer_command_from_text(text)
    if inferred_dict is None:
        return None

    prefix = inferred_dict['prefix']
    suffix = inferred_dict['suffix']
    category = inferred_dict['category']
    species = validate_correct_species(inferred_dict['species'])
    new_dict = {
        'prefix': prefix,
        'suffix': suffix,
        'category': category,
        'species': species,
    }

    suffix_type = phrase_dict[category]['suffix_type']
    if suffix_type == 'genestring':
        suffix_corrected = validate_correct_genestr(suffix, species=species)
        question = 'gene_string'
    elif suffix_type == 'regionsinglestring':
        suffix_corrected = validate_correct_single_region(suffix, species=species)
        question = 'regionsingle_string'
    elif suffix_type == 'regionmultistring':
        suffix_corrected = validate_correct_multiregion(suffix, species=species)
        question = 'regionmulti_string'
    elif suffix_type == 'celltypestring':
        suffix_corrected = validate_correct_celltypestr(suffix)
        question = 'celltype_string'
    elif suffix_type == 'celltype_dataset_timepoint_string':
        suffix_corrected = validate_correct_celltypedatasettimepoint(
                suffix)
        question = 'celltype_dataset_timepoint_string'
    elif suffix_type == 'timepoint':
        suffix_corrected = validate_correct_timepoint(suffix)
        question = 'timepoint'
    elif suffix_type == "species_genestring":
        suffix_corrected = validate_correct_species_genestr(species, suffix)
        question = 'species_gene_string'
    else:
        raise ValueError('Category not implemented')

    new_dict.update({
        'suffix_corrected': suffix_corrected,
        'question': question,
    })
    return new_dict
