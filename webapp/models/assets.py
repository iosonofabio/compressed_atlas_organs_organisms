from config import configuration as config

feature_types = config['feature_types']
fn_atlasd = config['paths']['compressed_atlas']
fn_GO = config['paths']['pathways']['GO']['mouse']

pseudocountd = {
    'gene_expression': 0.5,
    'chromatin_accessibility': 0.01,
}
