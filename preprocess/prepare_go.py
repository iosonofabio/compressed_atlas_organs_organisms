# vim: fdm=indent
'''
author:     Fabio Zanini
date:       24/05/22
content:    Prepare gene ontology (GO) categories for querying.
'''
import os
import sys
import gzip
import h5py
import pandas as pd


data_fdn = '../webapp/static/scData/'


if __name__ == '__main__':

    print('Connecting GO codes to gene lists')
    fn_go_genes = '../data/gene_ontology/mouse_mgi.gaf.gz'
    df = pd.read_csv(fn_go_genes, sep='\t', skiprows=36, header=None, usecols=[2, 3, 4])
    df.columns = ['Gene', 'relation', 'GO']

    df = df.loc[df['relation'] == 'involved_in']

    print('Fwd dic')
    dic = {}
    for key, val in df.groupby('GO'):
        tmp = list(set(list(val['Gene'].fillna('').values)) - set(''))
        if tmp:
            dic[key] = tmp

    print('Reverse dic')
    dic_inv = {}
    for key, val in df.groupby('Gene'):
        tmp = list(set(list(val['GO'].fillna('').values)) - set(''))
        if tmp:
            dic_inv[key] = tmp

    print('Giving names to GO codes')
    fn_go_obo = '../data/gene_ontology/go.obo.gz'
    dic_names = {}
    dic_ids = {}
    with gzip.open(fn_go_obo, 'rt') as f:
        goid = None
        goname = None
        for line in f:
            if line == '[Term]\n':
                goid = None
                goname = None
                continue
            elif line == '\n':
                if goname is not None:
                    dic_names[goname] = goid
                    dic_ids[goid] = goname
                continue
            elif line.startswith('id: GO:'):
                goid = line[4:-1]
                continue
            elif line.startswith('name: '):
                goname = line[6:-1]
                continue

    print('Connect the two')
    goids = list(dic.keys())
    for goid in goids:
        if goid in dic_ids:
            dic[dic_ids[goid]] = dic[goid]

    # Generic GO term, skip them
    goname_blacklist = (
        'biological_process',
    )
    dic_inv_names = {}
    for gene, goids in dic_inv.items():
        tmp = []
        for goid in goids:
            if goid in dic_ids:
                goname = dic_ids[goid]
                if goname in goname_blacklist:
                    continue
                tmp.append(goname)

        # skip genes with zero non-blacklisted GO terms
        if tmp:
            dic_inv_names[gene] = tmp

    print('Store to file (h5 not ideal here because of the jagged arrays)')
    import pickle
    fn_out = data_fdn + 'mouse_GO_tables.pkl'
    with open(fn_out, 'wb') as fout:
        pickle.dump({
            'GO_to_genes': dic,
            'genes_to_GO': dic_inv_names,
        }, fout)
