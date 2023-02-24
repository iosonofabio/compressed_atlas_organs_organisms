import pickle
import h5py
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from urllib.parse import urlencode

from config import configuration as config
from validation.celltypes import (
        validate_correct_celltypestr,
        )
from validation.differential_measurement import get_differential_conditions
from models.lazy import (
        get_tissues,
        get_celltypes,
        get_feature_orderd,
        get_feature_annotationd,
        get_feature_annotation,
        get_gene_annotationd,
        get_region_annotationd,
        get_gene_matrixd,
        )
from models.assets import (
        feature_types,
        fn_atlasd,
        fn_GO,
        pseudocountd,
        )



def get_counts(
        df_type,
        features=None,
        species='mouse',
        missing='throw',
        key="average",
        feature_type='gene_expression',
        ):
    # TODO: improve this docstring!!
    '''Read the h5 file with the compressed atlas


    gives the name of dataset we want as an input
    celltype / celltype_dataset / celltype_dataset_timepoint

    Args:
        key: Whether to return the gene expression average or the fraction
          of expressing cells.
    '''
    feature_order = get_feature_orderd(feature_type, species)

    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        columns = h5_data[feature_type][df_type]['index'].asstr()[:]

        # Lazy structure, for speed
        counts = h5_data[feature_type][df_type][key]
        if features is not None:
            index = []
            for gene in features:
                # Try perfect match
                if gene in feature_order.index:
                    if gene not in index:
                        index.append(gene)
                    continue

                # Try imperfect match
                candidates = feature_order.index[feature_order.index.str.match(gene)]
                if len(candidates) == 0:
                    if missing == 'throw':
                        raise KeyError('Gene match not found: '+gene)
                    else:
                        continue
                for cand in candidates:
                    if cand not in index:
                        index.append(cand)

            idx = feature_order.loc[index].values
            # Slices of h5 must be ordered...
            tmp = pd.Series(idx, index=np.arange(len(idx))).to_frame(name='idx')
            tmp = tmp.sort_values('idx')
            tmp['new_order'] = np.arange(len(tmp))
            counts = np.array(counts[:, tmp['idx'].values])
            # ... now reorder them according to the user's needs
            counts = counts[:, tmp.sort_index()['new_order'].values]
        else:
            index = feature_order.index
        # NOTE: Should be float32 already, check
        counts = counts[:].astype(np.float32)

    dataframe = pd.DataFrame(
            data=counts.T,
            index=index,
            columns=columns,
            )

    return dataframe


def get_number_cells(
        df_type,
        species='mouse',
        feature_type='gene_expression',
    ):
    '''Read the number of cells per time point per dataset from file'''
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas) as f:
        dic = f[feature_type][df_type]
        labels = dic['index'].asstr()[:]
        values = dic['cell_count'][:]
        ncells = pd.Series(data=values, index=labels)

    return ncells


def get_data_overtime_1feature(
        feature,
        species='mouse',
        feature_type='gene_expression',
    ):
    '''Get JS plotly code for big heatmap

    NOTE: this technically breaks the API concept. Let's keep it in mind
    and see what we can do. The positive is that a bunch of computations
    happen on the server... wait, is that good?
    '''
    from plotly import colors as pcolors

    ncells = get_number_cells(
            'celltype_dataset_timepoint',
            species=species,
            feature_type=feature_type,
            )
    countg = get_counts(
            'celltype_dataset_timepoint', features=[feature],
            species=species,
            feature_type=feature_type,
            ).iloc[0]

    # FIXME: this should not be necessary!
    countg = countg.fillna(0)

    # Sort the rows
    timepoint_order = config['order']['timepoint'][species]
    dataset_order = config['order']['dataset'][species]

    timepoint_orderd = {x: i for i, x in enumerate(timepoint_order)}
    dataset_orderd = {x: i for i, x in enumerate(dataset_order)}
    ncells.index = ncells.index.str.split('_', 1, expand=True)
    ncells = ncells.unstack(0, fill_value=0)
    ncells[['dataset', 'timepoint']] = ncells.index.str.split('_', expand=True).tolist()
    ncells['timepoint_order'] = ncells['timepoint'].map(timepoint_orderd)
    ncells['dataset_order'] = ncells['dataset'].map(dataset_orderd)
    ncells.sort_values(['timepoint_order', 'dataset_order'], inplace=True)
    countg.index = countg.index.str.split('_', 1, expand=True)
    countg = countg.unstack(0, fill_value=-1).loc[ncells.index]
    ncells = ncells.loc[:, countg.columns]

    # Sort the columns
    distance = pdist(countg.T.values)
    Z = linkage(distance, optimal_ordering=True)
    new_order = leaves_list(Z)

    celltypes = countg.columns.tolist()
    celltypes_hierarchical = countg.columns[new_order].tolist()
    row_labels = ncells.index.tolist()

    xticks = list(countg.columns)
    yticks = list(ncells.index)
    yticktext = []
    for i, label in enumerate(ncells.index):
        tp = label.split('_')[1]
        if (i == 0) or (tp != yticktext[-1]):
            yticktext.append(tp)
        else:
            yticktext.append('')

    return {
        'feature_type': feature_type,
        'feature': feature,
        'measurement': countg.T.to_dict(),
        'ncells': ncells.T.to_dict(),
        'celltypes': celltypes,
        'celltypes_hierarchical': celltypes_hierarchical,
        'row_labels': row_labels,
        'xticks': xticks,
        'yticks': yticks,
        'yticktext': yticktext,
        }


def get_data_overtime_1celltype(
        celltype,
        features,
        species='mouse',
        feature_type='gene_expression',
    ):
    '''Get JS plotly code for big heatmap

    NOTE: this technically breaks the API concept. Let's keep it in mind
    and see what we can do. The positive is that a bunch of computations
    happen on the server... wait, is that good?
    '''

    ncells = get_number_cells(
            'celltype_dataset_timepoint',
            species=species,
            feature_type=feature_type,
            )
    countg = get_counts(
            'celltype_dataset_timepoint', features=features,
            species=species,
            feature_type=feature_type,
            ).T

    # Only select rows with the correct cell type
    idx = []
    for label in countg.index:
        celltype_row = label.split('_')[0]
        celltype_row_validated = validate_correct_celltypestr(celltype_row)
        if celltype_row_validated == celltype:
            idx.append(label)
    countg = countg.loc[idx]

    # Same selection for the number of cells
    idx = []
    for label in ncells.index:
        celltype_row = label.split('_')[0]
        celltype_row_validated = validate_correct_celltypestr(celltype_row)
        if celltype_row_validated == celltype:
            idx.append(label)
    ncells = ncells.loc[idx]

    # Sort the rows
    timepoint_order = config['order']['timepoint'][species]
    dataset_order = config['order']['dataset'][species]

    countg['_n_cells'] = ncells.loc[countg.index]
    countg['_idx1'] = [timepoint_order.index(x.split('_')[2]) for x in countg.index]
    countg['_idx2'] = [dataset_order.index(x.split('_')[1]) for x in countg.index]
    countg = countg.sort_values(['_idx1', '_idx2'])
    ncells = countg['_n_cells']
    del countg['_idx1']
    del countg['_idx2']
    del countg['_n_cells']

    # Sort the columns
    if len(features) <= 2:
        new_order = list(range(len(features)))
    else:
        distance = pdist(countg.T.values)
        Z = linkage(distance, optimal_ordering=True)
        new_order = leaves_list(Z)

    features = countg.columns.tolist()
    features_hierarchical = countg.columns[new_order].tolist()
    row_labels = ncells.index.tolist()

    xticks = list(countg.columns)
    yticks = list(ncells.index)
    yticktext = []
    for i, label in enumerate(ncells.index):
        tp = label.split('_')[2]
        if (i == 0) or (tp != yticktext[-1]):
            yticktext.append(tp)
        else:
            yticktext.append('')

    return {
        'feature_type': feature_type,
        'celltype': celltype,
        'features': features,
        'features_hierarchical': features_hierarchical,
        'measurement': countg.T.to_dict(),
        'ncells': ncells.T.to_dict(),
        'row_labels': row_labels,
        'xticks': xticks,
        'yticks': yticks,
        'yticktext': yticktext,
        }


def get_features_correlated(
        features,
        species="mouse",
        feature_type='gene_expression',
        correlates_type='gene_expression',
        data_type='celltype',
        chain=True,
        correlation_min=0.15,
        nmax=None,
        ):
    '''Get features that are "friends" (correlated) with a list of features'''
    if len(features) == 0:
        if chain:
            return ''
        return {}

    # Compute on the fly, it's around 15ms per feature
    features_all = get_feature_orderd(feature_type, species)
    idx = features_all.loc[features].values

    counts_fea = get_counts(
        data_type,
        feature_type=feature_type,
        species=species,
    ).values.T
    countsc_fea = counts_fea - counts_fea.mean(axis=0)
    sigma_fea = counts_fea.std(axis=0)
    if feature_type == correlates_type:
        features_corr = features_all
        counts_corr = counts_fea
        countsc_corr = countsc_fea
        sigma_corr = sigma_fea
    else:
        features_corr = get_feature_orderd(correlates_type, species)
        counts_corr = get_counts(
            data_type,
            feature_type=correlates_type,
            species=species,
        ).values.T
        countsc_corr = counts_corr - counts_corr.mean(axis=0)
        sigma_corr = counts_corr.std(axis=0)

    if chain:
        correlates = []
    else:
        correlates = {}
    for i in idx:
        corr_i = (countsc_fea[:, i] * countsc_corr.T).mean(axis=1)
        den_i = (sigma_fea[i] * sigma_corr.T) + 1e-8
        corr_i /= den_i
        corr_i[np.isnan(corr_i)] = 0

        # Self is always the highest
        if feature_type == correlates_type:
            itop = np.argsort(corr_i)[::-1][1:6]
        else:
            itop = np.argsort(corr_i)[::-1][:5]

        # Cutoff
        tmp = corr_i[itop]
        tmp = tmp[tmp >= correlation_min]
        itop = itop[:len(tmp)]

        correlates_i = features_corr.index[itop]

        if nmax is not None:
            correlates_i = correlates_i[:nmax]

        if chain:
            for feature in correlates_i:
                if feature not in correlates:
                    correlates.append(feature)

        else:
            correlates[features[i]] = np.array(correlates_i)

    if chain:
        return ",".join(correlates)
    return correlates


def get_features_nearby(
        features,
        species="mouse",
        target_types=['gene_expression'],
        chain=True,
        dmax=50000,
        include_query=False,
        ):
    '''Get features that are "friends" (correlated) with a list of features'''
    if isinstance(features, str):
        features = features.split(',')

    if len(features) == 0:
        if chain:
            return ''
        return {}

    if chain:
        targets = []
    else:
        targets = {x: [] for x in features}

    for target_type in target_types:
        target_annos = get_feature_annotationd(target_type, species)
        chroms = target_annos['chrom']

        # NOTE: Poor man's fix to inconsistencies in chromosome names
        if chroms.str.startswith('chr').sum() > 0:
            chr_hack = 'add'
        else:
            chr_hack = 'remove'

        for feature in features:
            feature_annod = get_feature_annotation(feature, species)
            feature_type = feature_annod['feature_type']
            feature_anno = feature_annod['annotation']
            chrom = feature_anno['chrom']

            if (chr_hack == 'add') and (not chrom.startswith('chr')):
                chrom = 'chr' + chrom
            elif (chr_hack == 'remove') and chrom.startswith('chr'):
                chrom = chrom[3:]

            target_annos_chrom = target_annos.loc[chroms == chrom]

            if len(target_annos_chrom) == 0:
                targets_i = []
            else:
                start = feature_anno['start']
                end = feature_anno['end']
                starts = target_annos_chrom['start']
                ends = target_annos_chrom['end']
                overlap_i = (starts < end) & (ends > start)

                overlap_i = target_annos_chrom.index[overlap_i].tolist()

                dist_i = np.min([
                        (starts - end).abs(),
                        (ends - start).abs(),
                        ], axis=0)

                # Self is always the highest
                if feature_type == target_type:
                    itop = np.argsort(dist_i)[1:6]
                else:
                    itop = np.argsort(dist_i)[:5]

                # Cutoff
                tmp = dist_i[itop]
                tmp = tmp[tmp <= dmax]
                itop = itop[:len(tmp)]

                targets_i = overlap_i + list(target_annos_chrom.index[itop])

            if chain:
                for feature in targets_i:
                    if feature not in targets:
                        targets.append(feature)
            else:
                targets[feature].extend(targets_i)

    if include_query:
        if chain:
            targets = features + targets
        else:
            for feature in features:
                targets[feature] = np.array([feature] + list(targets[feature]))

    if chain:
        return ",".join(targets)
    return targets


def get_marker_features(
        celltypes,
        species='mouse',
        feature_type='gene_expression',
        chain=True,
        ntop=10):
    '''Get markers for cell types'''

    celltype_dict_inv = {val: key for key, val in celltype_dict.items()}
    celltype_dict_inv.update({ct: ct for ct in celltype_dict_inv.values()})

    # Sometimes we get a string
    if isinstance(celltypes, str):
        celltypes = celltypes.split(',')

    if feature_type == 'gene_expression':
        key = 'fraction'
    else:
        key = 'average'

    features_all = get_feature_orderd(feature_type, species)
    counts = get_counts(
        'celltype',
        feature_type=feature_type,
        species=species,
        key=key,
    )
    if chain:
        markers = []
    else:
        markers = {}

    # Go over the requested cell types
    for celltype in celltypes:
        if celltype not in counts.columns:
            return None

        val_ct = counts[celltype].values

        other_celltypes = [x for x in counts.columns if x != celltype]
        val_other = counts[other_celltypes].values

        diff_closest = val_ct - val_other.max(axis=1)
        diff_avg = val_ct - val_other.mean(axis=1)

        # Candidate markers: either up compared to closest, or compared to
        # average cell type
        idx_cand = list(
            set(np.argsort(diff_closest)[-ntop:]) | set(np.argsort(diff_avg)[-ntop:]),
        )
        cands = list(features_all.iloc[idx_cand].index)

        if chain:
            for cand in cands:
                if cand not in markers:
                    markers.append(cand)
        else:
            markers[celltype] = cands

    if chain:
        return ",".join(markers)
    return markers


def get_differential_url(
        suffix,
        feature_type='gene_expression',
        kind='both',
        n_features=20):
    '''Get differential expression url for up-, downregulated features, or both'''
    baseline = config.get('baseline', 'normal').lower()

    diff_dict = get_differential_conditions(suffix)
    diff_request = {
        'feature_type': feature_type,
        'comparison': diff_dict['kind'],
        'ct1': diff_dict['conditions'][0]['celltype'],
        'ct2': diff_dict['conditions'][1]['celltype'],
        'ds1': diff_dict['conditions'][0]['dataset'],
        'ds2': diff_dict['conditions'][1]['dataset'],
        'tp1': diff_dict['conditions'][0]['timepoint'],
        'tp2': diff_dict['conditions'][1]['timepoint'],
        'dis1': diff_dict['conditions'][0].get('condition', baseline),
        'dis2': diff_dict['conditions'][1].get('condition', baseline),
        'kind': kind,
        'n_features': n_features,
    }
    url = urlencode(diff_request)
    return url


def get_data_differential(
        conditions,
        kind='both',
        n_features=20,
        features=None,
        feature_type='gene_expression',
        species='mouse',
        ):
    '''Get differentially expressed, up- or downregulated features

    the special phease " in " separated the split (e.g. disease) from
    the celltype/timepoint/dataset specification (e.g. basophil ACZ P7)
    '''
    condition_name = config['condition'].lower()
    baseline_name = config.get('baseline', 'normal').lower()

    # If requested, compute DEGs on the fly
    if features is None:
        dfs = []

    # In any case, collect data to show in the heatmap
    dfs_alltypes = []
    for condition in conditions:
        if condition.get('condition', baseline_name) == condition_name:
            dfi = get_counts(
             'celltype_dataset_timepoint_disease',
             features=features,
             feature_type=feature_type,
             species=species,
             )
        else:
            dfi = get_counts(
                'celltype_dataset_timepoint',
                features=features,
                feature_type=feature_type,
                species=species,
                )

        # Get data for only this condition if computation of DEGs is requested
        if features is None:
            idx_col = '_'.join([
                condition['celltype'],
                condition['dataset'],
                condition['timepoint'],
                ])
            if idx_col not in dfi.columns:
                raise ValueError(f'Condition not found: {idx_col}')
            dfi_cond = dfi.loc[:, idx_col]
            dfs.append(dfi_cond)

        # In any case, get data for this condition, but all cell types
        idx_col = dfi.columns.str.contains(
                '_'+condition['dataset']+'_'+condition['timepoint'])
        dfi_alltypes = dfi.loc[:, idx_col]

        # Rename the columns as the cell types
        dfi_alltypes.columns = [x.split('_')[0] for x in dfi_alltypes.columns]
        dfs_alltypes.append(dfi_alltypes)

    if features is None:
        # Get log2 fold change
        log2fc = np.log2(dfs[0] + pseudocountd[feature_type]) - \
                np.log2(dfs[1] + pseudocountd[feature_type])

        # Get DEGs
        dmg_up = log2fc.nlargest(n_features).index.tolist()
        dmg_down = log2fc.nsmallest(n_features).index.tolist()

        if kind == 'up':
            features = dmg_up
        elif kind == 'down':
            features = dmg_down
        else:
            features = dmg_up
            for fea in dmg_down[::-1]:
                if fea not in features:
                    features.append(fea)

    # Use or recycle data to make heatmap
    # FIXME: make an expected total ordering including disease
    celltypes = dfs_alltypes[0].columns.tolist()
    for celltype in dfs_alltypes[1].columns:
        if celltype not in celltypes:
            celltypes.append(celltype)
    dfs_alltypes = [dfi.loc[features] for dfi in dfs_alltypes]
    for i, dfi in enumerate(dfs_alltypes):
        for celltype in celltypes:
            if celltype not in dfi.columns:
                dfi[celltype] = 0
        dfs_alltypes[i] = dfi

    return dfs_alltypes


def get_data_condition(
        features=None,
        feature_type='gene_expression',
        species='mouse',
        datasets=None,
        ):
    '''Get heatmap data for non-baseline conditions, e.g. disease

    You can limit it to specific datasets with the keyword arg.
    '''
    df_ho = get_counts(
        'celltype_dataset_timepoint_disease',
        features=features,
        feature_type=feature_type,
        species=species,
        )

    df_normal = get_counts(
            'celltype_dataset_timepoint',
            features=features,
            feature_type=feature_type,
            species=species,
        )
    # Restrict to disease celltypes, datasets, and timepoints
    # NOTE: no dataset took disease from a timepoint that has no normal.
    # However, some cell types might be condition-specific
    for key in df_ho:
        if key not in df_normal.columns:
            # Default to zero expression
            df_normal[key] = 0
    df_normal = df_normal.loc[df_ho.index, df_ho.columns]

    # Split by dataset and timepoint, and unstack
    df_dict = {'data': df_ho.T, 'data_baseline': df_normal.T}
    result = []

    # Check which datasets have disease
    dataset_timepoints = config['defaults']['condition']['dataset_timepoint']
    datasets_chosen = []
    timepointd = {}
    for dstp in dataset_timepoints:
        dataset, timepoint = dstp.split('_')
        # If restricted to certain datasets, skip the rest
        if (datasets is not None) and (dataset not in datasets):
            continue
        if dataset not in datasets_chosen:
            datasets_chosen.append(dataset)
        if dataset not in timepointd:
            timepointd[dataset] = []
        timepointd[dataset].append(timepoint)

    for dataset in datasets_chosen:
        for timepoint in timepointd[dataset]:
            item = {
                'dataset': dataset,
                'timepoint': timepoint,
            }
            for key, df in df_dict.items():
                dfi = df.loc[df.index.str.endswith(f'{dataset}_{timepoint}')]
                dfi.index = dfi.index.str.split('_', expand=True).get_level_values(0)
                item[key] = dfi
            result.append(item)

    return result


def get_celltype_abundances(
        timepoints=None,
        kind='qualitative',
        species='mouse'):
    '''Get cell type abundances at a certain time point'''
    ncells_tp = get_number_cells('celltype_dataset_timepoint')

    if timepoints is not None:
        idx = []
        for timepoint in timepoints:
            idx_tp = ncells_tp.index[ncells_tp.index.str.contains('_'+timepoint)]
            idx.extend(list(idx_tp))
        ncells_tp = ncells_tp.loc[idx]

    ncells_tp.index = [x.split('_')[0] for x in ncells_tp.index]

    ntot = ncells_tp.sum()
    fracs = 1.0 * ncells_tp / ntot

    if kind == 'quantitative':
        return fracs

    # Sort them naturally
    fracs.sort_values(ascending=False, inplace=True)

    result = {
        'major': [],
        'minor': [],
        'rare': [],
        'missing': [],
    }
    for ct, fr in fracs.items():
        if fr > 0.1:
            result['major'].append(ct)
        elif fr > 0.03:
            result['minor'].append(ct)
        elif fr > 0:
            result['rare'].append(ct)
        else:
            result['missing'].append(ct)

    return result


def get_gene_ids(genes, species='mouse'):
    '''Get the ids (MGI/GeneCards) of a list of genes'''
    if species == 'mouse':
        fn = config['paths']['gene_ids']
        df = pd.read_csv(fn, sep='\t', index_col=0)
        mgi_dict = df['MGI_id'].to_dict()
        human_dict = df['HumanGeneName'].fillna('').to_dict()
        id_dict = {}
        for gene in genes:
            new_name = human_dict.get(gene, '')
            if new_name == '':
                new_name = mgi_dict.get(gene, '')
            id_dict[gene] = new_name
    elif species == 'human':
        # GeneCards uses gene names as URL ids
        id_dict = {g: g for g in genes}
    elif species == 'lemur':
        # Seems like lemur and human share gene names
        return get_gene_ids(genes, species='human')
    else:
        raise NotImplementedError

    return id_dict


def get_gene_ontology_terms(genes, species='mouse'):
    '''Get the GO terms for these genes'''
    genesd = get_orthologs(genes, species, 'mouse')

    with open(fn_GO, 'rb') as fin:
        table = pickle.load(fin)['genes_to_GO']

    result = {}
    for gene_mouse, gene_orig in zip(genesd['mouse'], genesd[species]):
        if gene_mouse in table:
            result[gene_orig] = table[gene_mouse]

    return result


def get_genes_in_GO_term(go_term, species='mouse'):
    '''Get the genes in a GO term'''
    with open(fn_GO, 'rb') as fin:
        table = pickle.load(fin)['GO_to_genes']

    if go_term not in table:
        raise KeyError("GO term not found")

    genes = table[go_term]
    if species != 'mouse':
        genes = get_orthologs(genes, 'mouse', species)[species]

    return genes


def get_orthologs(genes, species, new_species):
    '''Connect orthologs from a species to another'''
    # Mouse lemur seems to use same gene names as human... is it a primate thing?
    if species == new_species:
        return {
            species: list(genes),
            new_species: list(genes),
            }
    elif species == 'lemur':
        dic = get_orthologs(genes, 'human', new_species)
        dic[species] = dic.pop('human')
        return dic
    elif new_species == 'lemur':
        dic = get_orthologs(genes, species, 'human')
        dic[new_species] = dic.pop('human')
        return dic
    elif (species == 'mouse') and (new_species == 'human'):
        fn = config['paths']['gene_ids']
        conv = pd.read_csv(fn, sep='\t', index_col=0, usecols=[0, 3])
        conv = conv.squeeze("columns")
        genes = [g for g in genes if g in conv.index]
        dic = conv.loc[genes].fillna('')
    elif (species == 'human') and (new_species == 'mouse'):
        fn = config['paths']['ortholog_genes']
        conv = pd.read_csv(fn, sep='\t', index_col=0)
        conv = conv.squeeze("columns")
        genes = [g for g in genes if g in conv.index]
        dic = conv.loc[genes].fillna('')
    else:
        raise NotImplementedError

    dic = dic[dic != '']
    return {
        species: list(dic.index),
        new_species: list(dic.values),
    }


def guess_genes_species(genes):
    '''Guess the species from the genes'''
    species = ('mouse',)
    for gene in genes:
        if len(gene) > 1:
            break
    else:
        return species

    if gene == gene.upper():
        species = ('human', 'lemur')
    return species


def get_data_species_comparison(species, species_baseline, genes):
    '''Get dataframe of average expression in two species.

    Genes with no ortholog will be ignored.
    '''
    gene_species = guess_genes_species(genes)

    # Several species could have the same
    if (species in gene_species) and (species_baseline in gene_species):
        genes_baseline = genes
    elif species in gene_species:
        dic = get_orthologs(genes, species, species_baseline)
        genes = dic[species]
        genes_baseline = dic[species_baseline]
    else:
        dic = get_orthologs(genes, species_baseline, species)
        genes = dic[species]
        genes_baseline = dic[species_baseline]

    # NOTE: Now genes are synched via orthology

    # Get both dataframes
    dfs = [
        get_counts('celltype', genes, species),
        get_counts('celltype', genes_baseline, species_baseline),
        ]

    # Restrict to common cell types
    celltypes_shared = set(dfs[0].columns.values) & set(dfs[1].columns.values)
    celltypes0 = [ct for ct in dfs[0].columns if ct in celltypes_shared]
    celltypes1 = [ct for ct in dfs[1].columns if ct in celltypes_shared]
    dfs = [dfs[0].loc[:, celltypes0], dfs[1].loc[:, celltypes1]]

    return dfs


def get_gsea(genes, species='mouse', gene_set='GO_Biological_Process_2021'):
    '''Get GSEA (gene set enrichment analysis)'''
    import gseapy as gp

    if gene_set == 'GO':
        gene_set = 'GO_Biological_Process_2021'
    elif gene_set == 'KEGG':
        if species == 'mouse':
            gene_set = 'KEGG_2019_Mouse'
        else:
            gene_set = 'KEGG_2021_Human'

    res = gp.enrichr(gene_list=genes,
                     gene_sets=[gene_set],
                     description='pathway',
                     outdir=None,
                     )

    res = res.results[['Term', 'Overlap', 'Adjusted P-value']]
    res.set_index('Term', inplace=True)

    return res


def get_kegg_urls(pathways):
    '''Get URLs to KEGG pathways by name'''
    kegg_df = pd.read_csv(
        config['paths']['pathways']['KEGG'], sep='\t', index_col='name')
    kegg_url = kegg_df['url'].to_dict()

    urls = [kegg_url.get(x, '') for x in pathways]
    return urls


def get_regions_related_to_gene(gene, relationships, species='mouse'):
    # TODO: Take a window ahead of TSS for now, we'll need to craft using the actual relationship

    regions = []

    if len(relationships) == 0:
        return regions

    if 'correlated' in relationships:
        regions_corr = get_features_correlated(
            [gene],
            feature_type='gene_expression',
            correlates_type='chromatin_accessibility',
        ).split(',')
        regions.extend(regions_corr)
        relationships.remove('correlated')

        if len(relationships) == 0:
            return regions

    gene_annotation = get_feature_annotationd(
            'gene_exression', species).loc[gene]
    chrom = gene_annotation['chrom']

    # Some genes/pseudogenes are not annotated by Ensembl
    if chrom == '':
        return regions

    region_annotations = get_feature_annotationd(
            'chromatin_accessibility', species)
    tss = gene_annotation['tss']
    strand = gene_annotation['strand']

    # FIXME: this is a hack - symptomatic for usage of a non-Ensembl genome
    # in the ATAC-Seq mapping. Need to chat with Emily about it because the
    # coordinates might be wrong
    chrom = 'chr' + chrom

    regions_chrom = region_annotations.loc[region_annotations['chrom'] == chrom]
    if len(regions_chrom) == 0:
        return regions

    # NB: gene/region start and end are always positive-strand
    gs = gene_annotation['start']
    ge = gene_annotation['end']
    rs = regions_chrom['start']
    re = regions_chrom['end']

    for rel in relationships:
        # Gene body
        if rel == 'gene body':
            idx = ((rs < gs) & (re > gs)) | ((rs >= gs) & (rs < ge))
        # Promoter (can be a single peak with exon 1)
        elif rel == 'promoter':
            if strand == 1:
                idx = (rs <= tss) & (re >= tss - 500)
            else:
                idx = (re >= tss) & (rs <= tss + 500)
        else:
            raise ValueError(f'Relationship not implemented: {rel}')

        regions_cand = regions_chrom.loc[idx].index.tolist()

        for region in regions_cand:
            if region not in regions:
                regions.append(region)
    return regions
