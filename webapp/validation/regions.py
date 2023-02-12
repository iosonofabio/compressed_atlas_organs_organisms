# vim: fdm=indent
'''
author:     Fabio Zanini
date:       21/12/22
content:    Validation of ATAC-Seq regions.
'''
import numpy as np
import pandas as pd

from models import (
    get_region_annotationd,
    get_regions_related_to_gene,
)
from validation.genes import validate_correct_gene



def validate_correct_single_region(regionstr, species='mouse'):
    '''Validate and correct misspellings for a single region name

    This excludes relationships with genes ATM since that might result in
    more than one region being selected.
    '''
    # TODO: allow for some mechanism to specify safe relationships

    # We generally want chrom_start_end
    parts = regionstr.split('-')

    if len(parts) != 3:
        return None

    # The second and this should be numbers
    chrom = parts[0]
    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError:
        return None

    if regionstr in get_region_annotationd(species).index:
        return regionstr

    # No perfect match, get closest match

    regions_matrixd = get_region_annotationd(species)
    regions_chrom = regions_matrixd.loc[regions_matrixd['chrom'] == chrom]

    # If not a known chromosome, nothing to do
    if len(regions_chrom) == 0:
        return None

    dstart = np.abs(regions_chrom['start'] - start)
    dend = np.abs(regions_chrom['end'] - end)
    ds_argmin = dstart.argmin()
    de_argmin = dend.argmin()
    if dstart[ds_argmin] < dend[de_argmin]:
        return regions_chrom.index[ds_argmin]
    else:
        return regions_chrom.index[de_argmin]


def validate_correct_multiregion(regionstr, species='mouse', missing='return_none'):
    '''Validate region names and correct misspellings if possible

    This function is quite complex because we want to leave the user the freedom
    to use intervals, e.g. chr1-1500k-1600k means "take all peaks on chromosome
    1, between coordinate 1,500,000 and 1,600,000". Morever several such
    intervals can be combined in a single query, e.g.

                    "chr1-1599k-8500k, chr3-8749k-12888k"

    Args:
        regionstr: see above.
        species: casual name for the species to look at: human, mouse, etc.
        missing: what to do with the missing features, e.g. if part of the
          query is not recognized. Accepted values are "return_none", which
          fails graciously; "throw", which fails with an exeption; and "ignore"
          which returns what was found and ignores the rest.

    Returns:
        If successful, validated string of regions matching the query. If failed,
        see above.
    '''

    # NOTE: This mess is mostly derived to accept voice recognition, not
    # touching for now.
    regions = regionstr.replace(' ', '').split(',')

    # Validate
    regionsv = []
    for region in regions:
        # We accept special syntax for gene-related queries, e.g. promoters, enhancers, suppressors
        if ':' in region:
            parts = region.split(':')
            # Other than 2 parts (relationship, gene) is unintelligible
            if len(parts) != 2:
                if missing == 'return_none':
                    return None
                elif missing == 'throw':
                    raise ValueError('Region not found:', region)
                else:
                    continue

            rel, gene = parts

            # Check relationship:
            # - p: (romoter)
            # - e: (nhancer)
            # - s: (uppressor)
            # - es: enhancer OR suppressor
            # - np: non-promoter, i.e. enhancers, suppressors, and other
            #       unclassified non-promoter regions
            # - gb: gene body, i.e. the transcript itself
            # - corr: peaks that are correlated with this gene
            if rel not in ('p', 'e', 's', 'es', 'np', 'gb', 'corr'):
                if missing == 'return_none':
                    return None
                elif missing == 'throw':
                    raise ValueError('Region not found:', region)
                else:
                    continue

            conv_table = {
                'p': ['promoter'], 'e': ['enhancer'], 's': ['suppressor'],
                'es': ['enhancer', 'suppressor'],
                'np': ['enhancer', 'suppressor', 'other_regulator'],
                'gb': ['gene body'],
                'corr': ['correlated'],
            }
            rels = conv_table[rel]

            # Now we can validate the actual gene name
            gene = validate_correct_gene(gene, species=species)

            # Ok, we actually need the region-gene relationships
            regions = get_regions_related_to_gene(gene, rels)

            if regions is None:
                if missing == 'return_none':
                    return None
                elif missing == 'throw':
                    raise ValueError('Region not found:', region)
                else:
                    continue
            
            # NOTE: we include the gene as well
            if gene not in regions:
                regionsv.append(gene)

            for region in regions:
                if region not in regionsv:
                    regionsv.append(region)

        # Else, it's a straight out region or an interval
        else:
            parts = region.split('-')

            # Other than 3 parts is unintelligible
            if len(parts) != 3:
                if missing == 'return_none':
                    return None
                elif missing == 'throw':
                    raise ValueError('Region not found:', region)
                else:
                    continue

            # Try exact match, we do this after splitting because it is expensive
            regions_matrixd = get_region_annotationd(species)
            if region in regions_matrixd.index:
                regionsv.append(region)
                continue

            # Inexact match, e.g. interval
            chrom, start_str, end_str = parts

            regions_chrom = regions_matrixd.loc[regions_matrixd['chrom'] == chrom]

            # If not a known chromosome, nothing to do
            if len(regions_chrom) == 0:
                if missing == 'return_none':
                    return None
                elif missing == 'throw':
                    raise ValueError('Region not found:', region)
                else:
                    continue

            # Flag for later
            failing = True

            coords_str = [start_str, end_str]
            coords = [-1, -1]
            for i, coord_str in enumerate(coords_str):
                # Exact match
                if coord_str.isdigit():
                    coords[i] = int(coord_str)
                    continue
                # Shorthand, e.g. 1900k
                if (len(coord_str) >= 2) and (coord_str[-1] in ('k', 'm')):
                    multiplier = 1000 if coord_str[-1] == 'k' else 1000000
                    try:
                        if '.' in coord_str[:-1]:
                            coords[i] = int(float(coord_str[:-1]) * multiplier)
                        else:
                            coords[i] = int(coord_str[:-1]) * multiplier
                        continue
                    except ValueError:
                        pass

                # Feature not found
                if missing == 'return_none':
                    return None
                elif missing == 'throw':
                    raise ValueError('Region not found:', region)
                else:
                    break
            else:
                failing = False

            if failing:
                continue

            # Take all regions in this interval
            start, end = coords

            idx = (regions_chrom['mid'] >= start) & (regions_chrom['mid'] <= end)
            regions_interval = regions_chrom.loc[idx]

            if len(regions_interval) == 0:
                if missing == 'return_none':
                    return None
                elif missing == 'throw':
                    raise ValueError('Region not found:', region)
                else:
                    continue

            # NOTE: this can get to be a *VERY* long list, so it'd be good to
            # keep the shorthand going down to the data models as slices. Too
            # messy given the current architecture, but not hard per se.
            for region in regions_interval.index:
                if region not in regionsv:
                    regionsv.append(region)

    regionstr = ','.join(regionsv)
    return regionstr

