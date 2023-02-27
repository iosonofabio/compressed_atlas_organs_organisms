# Web imports
from flask import request, jsonify, abort
from flask_restful import Resource

# Data import
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


# Helper functions
from config import configuration as config
from models import (
        get_speciess,
        get_counts,
        get_data_overtime_1feature,
        get_data_overtime_1celltype,
        get_data_condition,
        get_features_correlated,
        get_features_nearby,
        get_marker_features,
        get_celltype_abundances,
        get_data_differential,
        get_data_species_comparison,
        get_gene_ids,
        get_gene_ontology_terms,
        get_genes_in_GO_term,
        get_orthologs,
    )
from models.lazy import (
        get_feature_annotationd,
    )
from validation.genes import validate_correct_genestr
from validation.regions import (
        validate_correct_single_region,
        validate_correct_multiregion,
        )
from validation.feature_mixes import (
        validate_correct_feature_mix,
        validate_correct_single_feature,
        )
from validation.celltypes import validate_correct_celltypestr


class MeasurementByCelltype(Resource):
    '''API for the compressed atlas data, to be visualized as a heatmap'''
    def post(self):
        '''No data is actually posted, but this is uncached and can be longer

        In particular, when many genes are requested at once, we run into the
        length limitation of GET requests. The downside of using POST is that
        there's no caching, bookmarking, and such.
        '''
        args = request.form
        return self.get(args=args)

    def get(self, args=None):
        '''Get data to plot a table of measurements (e.g. GE) by cell type'''
        if args is None:
            args = request.args
        species = args.get("species")
        # Sometimes we are switching to a new species, if so get orthologs
        new_species = args.get("newSpecies")
        tissue = args.get("tissue")
        featurestring = args.get("feature_names")

        # A cap on gene names to avoid overload is reasonable
        featurestring = ','.join(featurestring.replace(' ', '').split(',')[:500])

        # Allow feature mixes, but keep them separate from now on
        # Each featurestring is not empty at the exit of the function below
        feature_stringd = validate_correct_feature_mix(
            featurestring, species=species,
        )
        if len(feature_stringd) == 0:
            print('no feature found!')
            return None

        result = {
            'data': [],
            'features': [],
            'feature_type': [],
            'features_hierarchical': [],
            'n_feature_types': len(feature_stringd),
            'data_fractions': [],
            'gene_ids': [],
            'feature_coords': [],
            'GO_terms': [],
        }

        # Store data to comupte hierarchical clustering of cell types
        dfl_for_ct_hierarchy = []

        for feature_type, featurestring in feature_stringd.items():
            # NOTE: this is where it gets tricky with canonical intervals
            feature_names = featurestring.split(',')

            # If we are switching species, get orthologs
            if new_species is not None:
                feature_names = get_orthologs(
                    feature_names, species, new_species,
                )[new_species]
                species = new_species
                missing_genes = 'skip'
            else:
                missing_genes = 'throw'

            try:
                df = get_counts(
                        "celltype",
                        feature_type=feature_type,
                        features=feature_names,
                        species=species,
                        tissue=tissue,
                        missing=missing_genes,
                        )
                if feature_type == 'gene_expression':
                    df_fractions = get_counts(
                            "celltype",
                            feature_type=feature_type,
                            features=feature_names,
                            key='fraction',
                            species=species,
                            tissue=tissue,
                            missing=missing_genes,
                            )

            except KeyError:
                print("Could not get counts from h5 file")
                return None

            # Just in case we skipped some
            feature_names = df.index.tolist()

            # Add gene ids and GO terms
            if feature_type == 'gene_expression':
                gene_ids = get_gene_ids(df.index, species=species)
                go_terms = get_gene_ontology_terms(feature_names, species=species)

            # Add feature coordinates
            feature_coords = get_feature_annotationd(
                    feature_type, species=species).loc[feature_names]
            feature_coordsstring = {}
            for feature_name in feature_names:
                chrom = feature_coords.at[feature_name, 'chrom']
                # FIXME: chr hack
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                start = feature_coords.at[feature_name, 'start']
                end = feature_coords.at[feature_name, 'end']
                fea_string = '-'.join([chrom, str(start), str(end)])
                feature_coordsstring[feature_name] = fea_string

            if feature_type == 'gene_expression':
                pseudocount = 0.5
            else:
                pseudocount = 0.01
            dfl = np.log10(df + pseudocount)

            if len(feature_names) <= 2:
                idx_features_hierarchical = list(range(len(feature_names)))
            else:
                idx_features_hierarchical = leaves_list(linkage(
                    pdist(dfl.values),
                    optimal_ordering=True),
                )
                idx_features_hierarchical = [int(x) for x in idx_features_hierarchical]

            result['data'].append(df.values.tolist())
            result['feature_type'].append(feature_type)
            result['features'].append(df.index.tolist())
            result['features_hierarchical'].append(idx_features_hierarchical)
            result['feature_coords'].append(feature_coordsstring)

            if feature_type == 'gene_expression':
                result['data_fractions'].append(df_fractions.values.tolist())
                result['gene_ids'].append(gene_ids)
                result['GO_terms'].append(go_terms)
            else:
                result['data_fractions'].append(None)
                result['gene_ids'].append(None)
                result['GO_terms'].append(None)

            dfl_for_ct_hierarchy.append(dfl.values)

        # Cell type hierarchical clustering is done on all features combined
        if len(dfl) == 1:
            dfl_for_ct_hierarchy = dfl_for_ct_hierarchy[0]
        else:
            dfl_for_ct_hierarchy = np.vstack(dfl_for_ct_hierarchy)

        idx_ct_hierarchical = leaves_list(linkage(
            pdist(dfl_for_ct_hierarchy.T),
            optimal_ordering=True),
        )
        idx_ct_hierarchical = [int(x) for x in idx_ct_hierarchical]

        result['celltypes'] = df.columns.tolist()
        result['celltypes_hierarchical'] = idx_ct_hierarchical
        result['species'] = species
        result['tissue'] = tissue

        return result


class MeasurementOvertime1Feature(Resource):
    '''Measurements (e.g. GE) for a single feature, across time'''
    def get(self):
        '''Get the measurement table (e.g. GE for plotting)'''
        species = request.args.get("species")
        new_species = request.args.get("newSpecies")
        tissue = request.args.get("tissue")

        feature_name = request.args.get("feature")
        try:
            featured = validate_correct_single_feature(
                feature_name, species=species,
            )
        except ValueError:
            return None
        feature_type = featured['feature_type']
        feature_name = featured['feature_name']

        # If we are switching species, get orthologs
        if new_species is not None:
            if feature_type == 'gene_expression':
                feature_name = get_orthologs(
                    [feature_name], species, new_species,
                )[new_species]

                if len(feature_name) == 0:
                    return None
            else:
                return None

            feature_name = feature_name[0]

            species = new_species
            missing_genes = 'skip'
        else:
            missing_genes = 'throw'

        data = get_data_overtime_1feature(
                feature_name,
                feature_type=feature_type,
                species=species,
                tissue=tissue,
                )

        if feature_type == "gene_expression":
            gene_id = get_gene_ids(
                    [feature_name], species=species)[feature_name]
            data['gene_id'] = gene_id

        data['similar_features'] = {}
        for correlates_type in config['feature_types']:
            similar_features = get_features_correlated(
                    [feature_name],
                    feature_type=feature_type,
                    correlates_type=correlates_type,
                    species=species,
                ).split(',')
            if (len(similar_features) > 0) and (similar_features[0] == feature_name):
                del similar_features[0]
            data['similar_features'][correlates_type] = similar_features

        return data


class MeasurementOvertime1Celltype(Resource):
    def get(self):
        species = request.args.get("species")

        celltype = request.args.get("celltype")
        celltype_validated = validate_correct_celltypestr(celltype)

        featurestring = request.args.get("feature_names")

        # A cap on gene names to avoid overload is reasonable
        featurestring = ','.join(featurestring.replace(' ', '').split(',')[:500])

        featured = validate_correct_feature_mix(
                featurestring,
                species=species)

        if len(featured) == 0:
            return None

        # If we are switching species, get orthologs
        new_species = request.args.get("newSpecies")
        if new_species is not None:
            # Converting chromatin features is impossible ATM
            if ('gene_expression' not in featured) or ('chromatin_accessibility' in featured):
                return None

            featured['gene_expression'] = ','.join(get_orthologs(
                    featured['gene_expression'].split(','), species, new_species,
                )[new_species])

            if len(featured['gene_expression']) == 0:
                return None

            species = new_species

        result = {
            'celltype': celltype,
            'ncells': [],
            'data': [],
            'features': [],
            'feature_type': [],
            'features_hierarchical': [],
            'n_feature_types': len(featured),
            'data_fractions': [],
            'gene_ids': [],
            'GO_terms': [],
        }

        for feature_type, featurestring in featured.items():
            feature_names = featurestring.split(',')
            data1t = get_data_overtime_1celltype(
                celltype_validated,
                feature_names,
                species=species,
                feature_type=feature_type,
            )

            result['features'].append(data1t['features'])
            result['features_hierarchical'].append(data1t['features_hierarchical'])
            result['data'].append(data1t['measurement'])
            result['ncells'].append(data1t['ncells'])
            # FIXME: this is a little tricky
            result['row_labels'] = data1t['row_labels']

            if feature_type == 'gene_expression':
                gene_ids = get_gene_ids(data1t['features'], species=species)
                go_terms = get_gene_ontology_terms(data1t['features'], species=species)
                result['gene_ids'].append(gene_ids)
                result['GO_terms'].append(go_terms)
            else:
                result['gene_ids'].append(None)
                result['GO_terms'].append(None)

        return result


class MeasurementCondition(Resource):
    '''API for condition vs normal data'''
    def get(self):
        '''Get measurement in a non-baseline condition, e.g. disease

        NOTE: the result from this API call differs from most others because
        this call can blend distinct data sets. Therefore, the returned data
        structure is a list of results, each of which describes a single
        combination of dataset, timepoint, and feature_type.

        The receiver (e.g. frontend) must take care to organise this data
        properly for plotting purposes.
        '''
        species = request.args.get("species")

        featurestring = request.args.get("feature_names")
        featured = validate_correct_feature_mix(
                featurestring,
                species=species)

        if len(featured) == 0:
            return None

        # Each feature type might have distinct time points
        # For the frontend menu, we want the union of all timepoints. When
        # clicked, only the relevant feature types are plotted
        result = []
        for feature_type, featurestring in featured.items():
            feature_names = featurestring.split(',')
            data1t = get_data_condition(
                features=feature_names,
                feature_type=feature_type,
                species=species,
            )

            # Each item is one dataset, timepoint, and feature type, and
            # includes both condition and baseline data
            for item in data1t:
                df = item['data']
                item['features'] = df.columns.tolist()
                item['celltypes'] = df.index.tolist()
                item['feature_type'] = feature_type

                if feature_type == 'gene_expression':
                    gene_ids = get_gene_ids(item['features'], species=species)
                    go_terms = get_gene_ontology_terms(
                            item['features'], species=species)
                    item['gene_ids'] = gene_ids
                    item['GO_terms'] = go_terms

                df_delta = np.log2(df + 0.5) - np.log2(item['data_baseline'] + 0.5)
                if len(item['celltypes']) > 1:
                    new_order = leaves_list(linkage(
                        pdist(df_delta.values),
                        optimal_ordering=True,
                        ))
                else:
                    new_order = [0]
                item['celltypes_hierarchical'] = df.index[new_order].tolist()
                if len(feature_names) > 1:
                    new_order = leaves_list(linkage(
                        pdist(df_delta.values.T),
                        optimal_ordering=True,
                        ))
                else:
                    new_order = [0]
                item['features_hierarchical'] = df.columns[new_order].tolist()

                item['data'] = item['data'].to_dict()
                item['data_baseline'] = item['data_baseline'].to_dict()

                result.append(item)

        return jsonify(result)


class MeasurementDifferential(Resource):
    '''API for generic differential expression'''
    def get(self):
        '''Get differential measurement data for plotting'''
        rqd = dict(request.args)
        conditions = [
                {'celltype': rqd['ct1'], 'timepoint': rqd['tp1'],
                 'dataset': rqd['ds1'], 'disease': rqd['dis1']},
                {'celltype': rqd['ct2'], 'timepoint': rqd['tp2'],
                 'dataset': rqd['ds2'], 'disease': rqd['dis2']},
        ]
        genes = rqd['genestr'].split(',')

        try:
            dfs = get_data_differential(conditions, genes=genes)
        except (KeyError, ValueError):
            return None

        # Get hierarchical clustering of cell types and genes
        df = dfs[0] - dfs[1]
        if len(genes) > 1:
            new_order = leaves_list(linkage(
                        pdist(df.values),
                        optimal_ordering=True,
                        ))
        else:
            new_order = [0]
        genes_hierarchical = df.index[new_order].tolist()
        new_order = leaves_list(linkage(
                    pdist(df.values.T),
                    optimal_ordering=True,
                    ))
        celltypes_hierarchical = df.columns[new_order].tolist()

        # Gene hyperlinks
        gene_ids = get_gene_ids(df.index)

        # Inject dfs into template
        heatmap_data = {
            'comparison': rqd['comparison'],
            'data': dfs[0].T.to_dict(),
            'data_baseline': dfs[1].T.to_dict(),
            'celltype': conditions[0]['celltype'],
            'celltype_baseline': conditions[1]['celltype'],
            'dataset': conditions[0]['dataset'],
            'dataset_baseline': conditions[1]['dataset'],
            'timepoint': conditions[0]['timepoint'],
            'timepoint_baseline': conditions[1]['timepoint'],
            'disease': conditions[0]['disease'],
            'disease_baseline': conditions[1]['disease'],
            'genes': dfs[0].index.tolist(),
            'celltypes': dfs[0].columns.tolist(),
            'genes_hierarchical': genes_hierarchical,
            'celltypes_hierarchical': celltypes_hierarchical,
            'gene_ids': gene_ids,
        }
        return heatmap_data


class MeasurementSpeciesComparison1Feature(Resource):
    '''Comparison between species'''
    def get(self):
        tissue = request.args.get('tissue')
        feature = request.args.get('feature').split(',')
        species_orig = request.args.get('species')
        feature_type = 'gene_expression'


        # Get the counts
        # NOTE: this function restricts to the intersection of cell types,
        # which makes the hierarchical clustering easy. In summary, both
        # genes and cell types are fully synched now
        speciess = get_speciess()

        dfd = {}
        for species in speciess:
            if species == species_orig:
                feature_orth = feature
            else:
                feature_orth = get_orthologs(
                    [feature], species_orig, species,
                )[species][0]

            df = get_counts(
                    "celltype",
                    feature_type=feature_type,
                    features=[feature_orth],
                    species=species,
                    tissue=tissue,
                    missing='throw',
                    )

            print(df)

            dfd[species] = df

        return {
            'data': {},  # TODO
            'celltypes': [],
            'celltypes_hierarchical': [],
            'speciess': list(speciess),
            'speciess_hierarchical': [],
            'feature_type': feature_type,
            'tissue': tissue,
        }

        ## Get hierarchical clustering of cell types
        ## NOTE: both dfs have the same celltypes (see above note)
        #new_order = leaves_list(linkage(
        #            pdist(df.values.T),
        #            optimal_ordering=True,
        #            ))
        #celltypes_hierarchical = df.columns[new_order].tolist()

        ## Gene hyperlinks (they hold for both)
        #gene_ids = get_gene_ids(df.index, species)

        ## Inject dfs into template
        ## NOTE: the whole converting DataFrame to dict of dict makes this quite
        ## a bit more heavy than it should be... just use a list of lists and
        ## accompanying lists of indices
        #result = {
        #    'data': dfs[0].T.to_dict(),
        #    'data_baseline': dfs[1].T.to_dict(),
        #    'genes': dfs[0].index.tolist(),
        #    'genes_baseline': dfs[1].index.tolist(),
        #    'celltypes': dfs[0].columns.tolist(),
        #    'celltypes_baseline': dfs[1].columns.tolist(),
        #    'genes_hierarchical': genes_hierarchical,
        #    'celltypes_hierarchical': celltypes_hierarchical,
        #    'genes_hierarchical_baseline': genes_hierarchical_baseline,
        #    'gene_ids': gene_ids,
        #    'species': species,
        #    'species_baseline': species_baseline,
        #}
        #return result


class PlotsForSeachGenes(Resource):
    def get(self):
        genestring = request.args.get("gene_names")
        gene_names = genestring.replace(' ', '').split(",")
        try:
            df = get_counts(
                    "celltype",
                    feature_type='gene_expression',
                    features=gene_names,
                ).T
        except KeyError:
            return None


        if len(gene_names) == 2:
            result = {}
            plot_df = df.filter(items=gene_names, axis=0)
            gene1 = plot_df.index[0]
            gene2 = plot_df.index[1]

            gene1_expr = list(plot_df.loc[gene1])
            gene2_expr = list(plot_df.loc[gene2])
            # plot_data = plot_df.to_json()
            result["gene1_name"] = gene1
            result["gene2_name"] = gene2
            result["gene1_expr"] = gene1_expr
            result["gene2_expr"] = gene2_expr
            result["cell_types"] = list(plot_df.columns)

        return result


class FeaturesCorrelated(Resource):
    def get(self):
        target_types = request.args.get("correlates_type")
        target_types = target_types.replace(" ", "").split(",")
        species = request.args.get("species")
        data_type = request.args.get("data_type", "celltype")
        n_target = int(request.args.get("n_correlates", 10))

        # Split features into feature types for convenience
        feature_names = request.args.get("feature_names")
        featured = validate_correct_feature_mix(
                feature_names,
                species=species,
        )

        # Collect all correlates of the requested types from all features
        features = []
        for feature_type, feature_names_type in featured.items():
            feature_names_type = feature_names_type.split(',')
            for target_type in target_types:
                featuresi = get_features_correlated(
                    feature_names_type,
                    species=species,
                    feature_type=feature_type,
                    correlates_type=target_type,
                    data_type=data_type,
                    chain=True,
                    nmax=n_target,
                    ).split(',')
                for feature in feature_names_type + featuresi:
                    if feature not in features:
                        features.append(feature)

        return ','.join(features)


class FeaturesNearby(Resource):
    def get(self):
        target_types = request.args.get("target_type")
        target_types = target_types.replace(" ", "").split(",")
        species = request.args.get("species")
        data_type = request.args.get("data_type", "celltype")
        distance_max = int(request.args.get("distance_max", 50000))

        # Split features into feature types for convenience
        feature_names = request.args.get("feature_names")
        featured = validate_correct_feature_mix(
                feature_names,
                species=species,
        )

        # Collect all correlates of the requested types from all features
        features = []
        for feature_type, feature_names_type in featured.items():
            featuresi = get_features_nearby(
                feature_names_type,
                species=species,
                target_types=target_types,
                chain=True,
                dmax=distance_max,
                include_query=True,
                ).split(',')
            for feature in featuresi:
                if feature not in features:
                    features.append(feature)

        return ','.join(features)


class GenesInGOTerm(Resource):
    def get(self):
        go_term = request.args.get("goTerm")
        species = request.args.get("species")
        genes = get_genes_in_GO_term(go_term, species)
        return ','.join(genes)


class CheckGenenames(Resource):
    def get(self):
        names = request.args.get("gene_names")
        names_validated = validate_correct_genestr(names)
        if names_validated is None:
            return {
                'outcome': 'fail',
                }
        else:
            return {
                'outcome': 'success',
                'genenames': names_validated,
                }


class MarkerGenes(Resource):
    def get(self):
        names = request.args.get("celltype_names")
        names_validated = validate_correct_celltypestr(names)
        if names_validated is None:
            return {
                'outcome': 'fail',
                }
        else:
            genenames = get_marker_features(names_validated)
            return {
                'outcome': 'success',
                'genenames': genenames,
                }


class CelltypeAbundance(Resource):
    def get(self):
        timepoint = request.args.get("timepoint")
        kind = request.args.get("kind")

        return {
            'outcome': 'success',
            'celltypeabundance': get_celltype_abundances(timepoint, kind=kind),
            }

