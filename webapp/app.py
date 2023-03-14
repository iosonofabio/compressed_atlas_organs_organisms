# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/05/22
content:    Main flask app for compressed atlases
'''
from flask import (
    Flask,
    send_from_directory,
    request,
    redirect,
    url_for,
    render_template as render_template_orig,
    abort,
)
from flask_restful import Api
from flask_cors import CORS
import numpy as np

from config import configuration as config
from api import (
    MeasurementByCelltype,
    MeasurementOvertime1Feature,
    MeasurementOvertime1Celltype,
    MeasurementCondition,
    PlotsForSeachGenes,
    FeaturesCorrelated,
    FeaturesNearby,
    GenesInGOTerm,
    MeasurementDifferential,
    MeasurementSpeciesComparison1Feature,
    MeasurementSpeciesComparison1Celltype,
    CheckGenenames,
    MarkerGenes,
    CelltypeAbundance,
)
from models import (
        get_speciess,
        get_tissues,
        get_celltypes,
        get_celltype_abundances,
        get_data_differential,
        get_gene_ids,
        get_features_correlated,
        get_data_species_comparison,
        get_orthologs,
        get_genes_in_GO_term,
        get_gsea,
        get_kegg_urls,
        pseudocountd,
)
from validation.genes import validate_correct_genestr
from validation.timepoints import validate_correct_timepoint
from text_recognition import mod as text_control_blueprint


##############################
# Contextualized templates
##############################
def render_template(*args, **kwargs):
    kwargs['species'] = kwargs.get('species', config['defaults']['species'])
    kwargs['speciess'] = get_speciess()
    kwargs['tissue'] = kwargs.get('tissue', config['defaults']['tissue'])
    kwargs['tissues'] = list(get_tissues(
            'gene_expression', kwargs['species']))
    kwargs['feature_types'] = config['feature_types']
    #kwargs['condition'] = config['condition'].capitalize()
    kwargs['unit_measurement'] = config['units']['counts']['long']
    kwargs['available_species'] = list(config['paths']['compressed_atlas'].keys())
    if 'short' in config['units']['counts']:
        kwargs['unit_measurement_short'] = config['units']['counts']['short']
    else:
        words = kwargs['unit_measurement'].split()
        kwargs['unit_measurement_short'] = ''.join([w[0] for w in words])
    return render_template_orig(
        *args, **kwargs,
    )


##############################
app = Flask(__name__, static_url_path="/static", template_folder="templates")
app_api = Api(app)
# NOTE: this might be unsafe
CORS(app)
with open('secret_key.txt') as f:
    app.config['SECRET_KEY'] = f.read()
##############################


##############################
# Views
##############################
@app.route("/")
def index():
    return redirect(url_for('text_control'))


# Control pages
@app.route("/start", methods=["GET"])
def text_control():
    """A single text bar to ask questions or post commands"""
    return render_template(
        "text_control.html",
        )


@app.route("/measurement_by_celltype", methods=['GET'])
def measurement_by_celltype():
    species = request.args.get('species')
    if species is None:
        species = config['defaults']['species']
    tissue = request.args.get('tissue')
    if tissue is None:
        tissue = config['defaults']['tissue']
    featurestring = request.args.get("featurestring")
    if featurestring is None:
        pathway = request.args.get("pathway")
        if pathway is not None:
            if ' (GO' in pathway:
                pathway = pathway[:pathway.find(' (GO')]
            try:
                features = get_genes_in_GO_term(pathway, species=species)
            except KeyError:
                return abort(404)

        else:
            features = config['defaults']['genestring'].split(',')
            if species in ('human', 'lemur'):
                features = get_orthologs(features, 'mouse', species)
        featurestring = ','.join(features)
    searchstring = featurestring.replace(" ", "")
    return render_template(
            "measurement_by_celltype.html",
            searchstring=searchstring,
            species=species,
            tissue=tissue,
            )


@app.route("/measurement_over_time_1feature", methods=['GET'])
def measurement_overtime_1feature():
    species = request.args.get('species')
    if species is None:
        species = config['defaults']['species']
    featurestring = request.args.get("featurestring")
    if featurestring is None:
        featurestring = config['defaults']['gene']
    searchstring = featurestring.replace(" ", "")

    return render_template(
            "measurement_overtime_1feature.html",
            searchstring=searchstring,
            species=species,
            )


@app.route("/measurement_over_time_1celltype", methods=['GET'])
def measurement_overtime_1celltype():
    species = request.args.get('species')
    if species is None:
        species = config['defaults']['species']
    tissue = request.args.get('tissue')
    if tissue is None:
        tissue = config['defaults']['tissue']
    celltype = request.args.get("celltype")
    if celltype is None:
        celltype = config['defaults']['celltype']
    genestring = request.args.get("genestring")
    if genestring is None:
        genestring = config['defaults']['genestring']
    searchstring = genestring.replace(" ", "")

    similar_genes = get_features_correlated(searchstring.split(',')).split(',')
    # Limit to a few
    similar_genes = similar_genes[:15]

    return render_template(
            "measurement_overtime_1celltype.html",
            celltype=celltype,
            searchstring=searchstring,
            species=species,
            tissue=tissue,
            similarGenes=similar_genes,
            # FIXME: this should depend on the species...
            celltypes=get_celltypes('gene_expression', species),
            )


@app.route("/measurement_condition", methods=["GET"])
def measurement_condition():
    """A sort of heatmap with condition (exercise vs not)"""
    species = request.args.get('species')
    if species is None:
        species = 'mouse'

    featurestring = request.args.get("featurestring")
    if featurestring is None:
        featurestring = config['defaults']['genestring']
    searchstring = featurestring.replace(" ", "")

    # Default dataset/timepoints combos
    dataset_timepoints = config['defaults']['condition']['dataset_timepoint']

    return render_template(
            "measurement_condition.html",
            searchstring=searchstring,
            datasetTimepoints=dataset_timepoints,
            species=species,
            )


# NOTE: This is where API and views break down and react would be better
@app.route("/heatmap_differential", methods=["GET"])
def heatmap_differential_measurement():
    """A sort of heatmap for differential measurment"""
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    species = request.args.get('species', 'mouse')

    rqd = dict(request.args)
    conditions = [
            {'celltype': rqd['ct1'],
             'timepoint': rqd['tp1'],
             'dataset': rqd['ds1'],
             'condition': rqd['dis1']},
            {'celltype': rqd['ct2'],
             'timepoint': rqd['tp2'],
             'dataset': rqd['ds2'],
             'condition': rqd['dis2']},
    ]

    # TODO: adapt to multiple feature types, i.e. feature_type="all"
    feature_type = rqd['feature_type']
    if feature_type == 'all':
        feature_types = config['feature_types']
    else:
         feature_types = [feature_type]

    plot_data = {
        'comparison': rqd['comparison'],
        'celltype': conditions[0]['celltype'],
        'celltype_baseline': conditions[1]['celltype'],
        'dataset': conditions[0]['dataset'],
        'dataset_baseline': conditions[1]['dataset'],
        'timepoint': conditions[0]['timepoint'],
        'timepoint_baseline': conditions[1]['timepoint'],
        'condition': conditions[0]['condition'],
        'condition_baseline': conditions[1]['condition'],
        'data': [],
        'data_baseline': [],
        'features': [],
        'features_hierarchical': [],
        'gene_ids': [],
        'feature_type': [],
    }

    for feature_type in feature_types:
        dfs = get_data_differential(
                conditions,
                kind=rqd['kind'],
                n_features=int(rqd['n_features']),
                feature_type=feature_type,
        )

        # Get hierarchical clustering of cell types and features
        # FIXME: this works poorly, trying out HC on the log
        df = np.log10(dfs[0] + pseudocountd[feature_type])
        new_order = leaves_list(linkage(
                    pdist(df.values),
                    optimal_ordering=True,
                    ))
        features_hierarchical = df.index[new_order].tolist()
        new_order = leaves_list(linkage(
                    pdist(df.values.T),
                    optimal_ordering=True,
                    ))
        celltypes_hierarchical = df.columns[new_order].tolist()

        # Gene hyperlinks
        if feature_type == 'gene_expression':
            gene_ids = get_gene_ids(df.index, species)
        else:
            gene_ids = []

        # NOTE: Inject dfs into template!! (why?)
        # FIXME: enable multiple feature types
        plot_data['data'].append(dfs[0].T.to_dict())
        plot_data['data_baseline'].append(dfs[1].T.to_dict())
        plot_data['features'].append(dfs[0].index.tolist())
        plot_data['celltypes'] = dfs[0].columns.tolist()
        plot_data['features_hierarchical'].append(features_hierarchical)
        plot_data['celltypes_hierarchical'] = celltypes_hierarchical
        plot_data['gene_ids'].append(gene_ids)
        plot_data['feature_type'].append(feature_type)

    # Set search string
    searchstring = ','.join(sum(plot_data['features'], []))

    return render_template(
            "measurement_differential.html",
            searchstring=searchstring,
            plotData=plot_data,
            )


@app.route("/list_celltypes/", methods=["GET"])
@app.route("/list_celltypes/<timepoint>", methods=["GET"])
def list_celltypes_timepoint(timepoint='', species="mouse"):
    '''List cell types and their abundances'''
    species = request.args.get('species', 'mouse')
    tissue = request.args.get(
            'tissue', get_tissues('gene_expression', species)[0])

    other_timepoints = ['Overall'] + config['order']['timepoint'][species]
    celltypes = list(get_celltypes(
        'gene_expression', species=species, tissue=tissue))

    if timepoint != '':
        timepoint = validate_correct_timepoint(timepoint)
        celltype_abundance = get_celltype_abundances(
                [timepoint],
                kind='quantitative',
                species=species,
                tissue=tissue,
                )
        celltypes_sorted = list(
                celltype_abundance.sort_values(ascending=False).index)
        celltype_abundance = celltype_abundance.to_dict()
        other_timepoints = [x for x in other_timepoints if x != timepoint]
    else:
        other_timepoints = [x for x in other_timepoints if x != 'Overall']
        celltypes_sorted = []
        celltype_abundance = ''

    return render_template(
            'list_celltypes.html',
            timepoint=timepoint,
            other_timepoints=other_timepoints,
            plotData={
                'celltype_abundance': celltype_abundance,
                'celltypes': celltypes,
                'celltypes_sorted': celltypes_sorted,
            },
            kind='quantitative',
            searchstring=timepoint,
            searchvalue="Timepoint",
            )


@app.route("/celltype_abundance/<timepoint>", methods=["GET"])
def plot_celltype_abundance(timepoint):
    '''Plot cell type abundances'''
    species = request.args.get('species', 'mouse')
    tissue = request.args.get(
            'tissue', get_tissues('gene_expression', species)[0])

    celltype_dict = get_celltype_abundances(
            timepoint,
            kind='quantitative',
            species=species,
            tissue=tissue,
            )
    return render_template(
            'celltype_abundance.html',
            timepoint=timepoint,
            celltypes=celltype_dict,
            searchstring=timepoint,
            )


@app.route("/barplot_gsea", methods=["GET", "POST"])
def plot_barplot_GSEA():
    '''Barplot for gene set enrichment analysis'''
    if request.method == "POST":
        args = request.form
    else:
        args = request.args

    genestring = args.get('genes')
    species = args.get('species')
    gene_set = args.get('gene_set', 'GO')

    genes = validate_correct_genestr(genestring, species=species).split(',')
    data = get_gsea(genes, species, gene_set=gene_set)

    if 'KEGG' in gene_set:
        pathway_urls = get_kegg_urls(data.index)
    else:
        pathway_urls = []

    # Cut too long results
    data = data.iloc[:15]

    return render_template(
        'barplot_gsea.html',
        species=species,
        plotData=dict(
            pathways=data.index.tolist(),
            pathways_urls=pathway_urls,
            overlap=data['Overlap'].values.tolist(),
            neglog10_p_value=(-np.log10(data['Adjusted P-value'].values)).tolist(),
            ),
        )


@app.route("/species_comparison/1feature", methods=["GET"])
def measurement_species_comparison_1feature():
    '''Compare measurement across multiple species, one feature'''
    feature = request.args.get('feature')
    if feature is None:
        feature = config['defaults']['gene']
    tissue = request.args.get('tissue')
    if tissue is None:
        tissue = config['defaults']['tissue']
    species_orig = request.args.get('species')
    if species_orig is None:
        species_orig = config['defaults']['species']

    searchstring = feature

    return render_template(
        'measurement_species_comparison_1feature.html',
        searchstring=searchstring,
        species=species_orig,
        tissue=tissue,
    )


@app.route("/species_comparison/1celltype", methods=["GET"])
def measurement_species_comparison_1celltype():
    '''Compare measurement across multiple species, one cell type'''
    celltype = request.args.get('celltype')
    if celltype is None:
        celltype = config['defaults']['celltype']
    tissue = request.args.get('tissue')
    if tissue is None:
        tissue = config['defaults']['tissue']
    species_orig = request.args.get('species')
    if species_orig is None:
        species_orig = config['defaults']['species']

    featurestring = request.args.get("featurestring")
    if featurestring is None:
        featurestring = config['defaults']['genestring']

    return render_template(
        'measurement_species_comparison_1celltype.html',
        searchstring=featurestring,
        species=species_orig,
        tissue=tissue,
        celltype=celltype,
    )



# Static assets (JS/CSS)
@app.route("/js/<path:path>")
def send_js(path):
    """JavaScript assets"""
    return send_from_directory("static/js", path)


@app.route("/css/<path:path>")
def send_css(path):
    """CSS stylesheets"""
    return send_from_directory("static/css", path)


@app.route("/favicon.ico")
def favicon():
    return send_from_directory("static", "favicon.ico")


# Blueprints
app.register_blueprint(text_control_blueprint)


# API endpoints
app_api.add_resource(MeasurementByCelltype, "/data/by_celltype")
app_api.add_resource(MeasurementOvertime1Feature, "/data/overtime_1feature")
app_api.add_resource(MeasurementOvertime1Celltype, "/data/overtime_1celltype")
app_api.add_resource(MeasurementCondition, "/data/condition")
app_api.add_resource(FeaturesCorrelated, "/data/features_correlated")
app_api.add_resource(FeaturesNearby, "/data/features_nearby")
app_api.add_resource(GenesInGOTerm, "/data/genes_in_go_term")
app_api.add_resource(MeasurementDifferential, "/data/differential")
app_api.add_resource(MeasurementSpeciesComparison1Feature,
                     "/data/speciescomparison/1feature")
app_api.add_resource(MeasurementSpeciesComparison1Celltype,
                     "/data/speciescomparison/1celltype")
app_api.add_resource(CheckGenenames, "/check_genenames")
app_api.add_resource(MarkerGenes, "/data/marker_genes")
app_api.add_resource(CelltypeAbundance, "/data/celltype_abundance")
# FIXME: this should not be a separate API endpoint
app_api.add_resource(PlotsForSeachGenes, "/data/by_celltype_2_genes")


# Main loop
if __name__ == "__main__":
    # app.run(debug=True)
    app.run(host="0.0.0.0", port=5000)
