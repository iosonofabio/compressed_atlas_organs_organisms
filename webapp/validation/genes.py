# vim: fdm=indent
'''
author:     Fabio Zanini
date:       30/04/22
content:    Validate gene names et al.
'''
from models import (
        get_feature_orderd,
        get_gene_matrixd
        )


def convert_numbers_in_gene_name(gene):
    '''Convert numbers in gene name into digits'''
    from text_recognition.assets import numbers

    # It's typical for gene names to end with a number (e.g. Car4)
    endswithdigit = False
    for i in range(len(gene)):
        tail = gene[len(gene) - 1 - i:]
        if tail.isdigit():
            endswithdigit = True
            continue
        break

    if endswithdigit:
        # No gene name is purely a number, this should be fine
        tail = tail[1:]
        ndigit = int(tail)
    else:
        # Check if we can convert the end to a digit
        for ntext, ndigit in numbers[::-1]:
            if gene.endswith(ntext):
                gene = gene[:-len(ntext)]+str(ndigit)
                endswithdigit = True
                break

    # If there are no convertible-to-digits at the end, we are done
    # FIXME: be more flexible than this, of course
    if not endswithdigit:
        return gene

    # If we found or converted an end digit, we should look for
    # internal digit-likes, e.g. Col1a1
    sfx = len(str(ndigit)) + 1
    for ntext, ndigit in numbers[::-1]:
        if gene[:-sfx].endswith(ntext):
            gene = gene[:-sfx-len(ntext)] + str(ndigit) + gene[-sfx:]
            break

    return gene


def validate_correct_gene(gene, species='mouse'):
    '''Validate and correct misspellings for a single gene name'''

    # Murine genes are capitalized
    # Human/monkey genes are upper
    if species == 'mouse':
        gene = gene.capitalize()
    else:
        gene = gene.upper()

    gene_order = get_feature_orderd('gene_expression', species)
    if gene in gene_order.index:
        return gene

    gene = convert_numbers_in_gene_name(gene)

    # Not found in whitelist, try to correct
    gene_matrix = get_gene_matrixd(species)
    gene_array = list(gene)[:gene_matrix.shape[1]]
    hamming = (gene_array != gene_matrix[:, :len(gene_array)]).sum(axis=1)

    # If the uncorrected gene is shorter (e.g. Col1) it can be a perfect match
    # for multiple whitelist genes (e.g. Col1a1, Col1a2), then ask for confirmation
    idx = (hamming == 0).nonzero()[0]
    # TODO: build mismatch scoring table based on natural English (e.g. p-t)
    if len(idx) == 0:
        idx = (hamming == 1).nonzero()[0]

    # If there is only one (partial) perfect match, take it. If multiple
    # perfect matches, ask. If no perfect and one imperfect match, take it.
    # If multiple imperfect or no imperfect, ask.
    if len(idx) == 1:
        gene_closest = gene_order.index[idx[0]]
        print('Getting closest gene:', gene, gene_closest)
        return gene_closest

    # At least one gene was not understood, ask for a written confirmation
    return None


def validate_correct_genestr(genestr, species='mouse', missing='return_none'):
    '''Validate gene names and correct misspellings if possible'''

    # TODO: check punctuation more accurately
    # Capitalization is taken care of in "validate_correct_gene"
    genes = genestr.strip(' ').replace('.', ',').replace(';', ',').split(',')

    # Validate
    genesv = []
    for gene in genes:
        genev = validate_correct_gene(gene, species=species)
        if genev is None:
            if missing == 'return_none':
                return None
            elif missing == 'throw':
                raise ValueError('Gene not found:', gene)
            else:
                continue
        genesv.append(genev)

    genestr = ','.join(genesv)
    return genestr
