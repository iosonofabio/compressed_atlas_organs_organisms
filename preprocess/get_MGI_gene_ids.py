# vim: fdm=indent
'''
author:     Fabio Zanini
date:       07/05/22
content:    Get dictionary connecting gene name and MGI gene id
'''
import os
import sys
import numpy as np
import pandas as pd


if __name__ == '__main__':

    fn = '../data/mouse_gene_names/MGI_Gene_Model_Coord.rpt'
    df = pd.read_csv(
            fn, sep='\t', usecols=[0, 1, 2],
            index_col=2, skiprows=1, header=None,
            )
    df.index.name = 'GeneName'
    df.columns = ['MGI_id', 'FeatureType']


    # Get human homolog and human gene name if present, because GeneCards
    # is better than MGI
    fn_homo = '../data/mouse_gene_names/HOM_MouseHumanSequence.rpt'
    df_homo = pd.read_csv(
            fn_homo, sep='\t',
            )
    dfi = df_homo[['DB Class Key', 'Common Organism Name', 'Symbol']].copy()
    orgname = dfi['Common Organism Name'].values
    dfi['Common Organism Name'] = ['mouse' if 'mouse' in x else 'human' for x in orgname]
    dfi.drop_duplicates(['DB Class Key', 'Common Organism Name'], inplace=True)
    dfi.set_index(['DB Class Key', 'Common Organism Name'], inplace=True)
    dfi = dfi['Symbol'].unstack(1, fill_value='')[['mouse', 'human']]
    dfid = dfi.set_index('mouse')['human'].to_dict()
    
    # mouse -> human also includes MGI, so it starts from df
    df['HumanGeneName'] = [dfid.get(x, '') for x in df.index]
    df.to_csv('../webapp/static/scData/mouse_gene_names.tsv', sep='\t', index=True)

    # human -> mouse has no MGI, so it takes dfi straight out  
    df_hm = dfi.drop_duplicates(subset='human').set_index('human')[['mouse']]
    df_hm.index.name = 'HumanGeneName'
    df_hm.columns = ['MouseGeneName']
    df_hm.to_csv(
            '../webapp/static/scData/human_mouse_gene_orthologs.tsv',
            sep='\t',
            index=True,
    )
