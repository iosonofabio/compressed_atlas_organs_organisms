# vim: fdm=indent
'''
author:     Fabio Zanini
date:       26/05/22
content:    Extract KEGG dict from scraping their website
'''
import numpy as np
import pandas as pd

data_fdn = '../webapp/static/scData/'


if __name__ == '__main__':

    # blacklists for weird organisms
    blacklist = ['fly', 'yeast', 'plant', 'HIV-1 New!']

    dic = {}
    fn = '../data/KEGG/kegg_pathways_scraping.txt'
    with open(fn) as f:
        line = 'initial'
        while line:
            line = f.readline()
            if line.startswith('0'):
                # Get id and type (useful for URL)
                pathway_id = line.rstrip('\n')
                ptype = ''
                if ' ' in pathway_id:
                    pathway_id, ptype = pathway_id.split(maxsplit=1)

                # Read description
                name = f.readline().rstrip('\n').strip()

                # Check blacklist
                for bit in blacklist:
                    if name.endswith(' - '+bit):
                        break
                else:
                    dic[pathway_id] = {'name': name, 'type': ptype}

    dic = pd.DataFrame(dic).T
    urls = []
    for pid, ptype in dic['type'].items():
        if 'M' in ptype:
            prefix = 'map'
        else:
            prefix = 'hsa'
        url = f'https://www.kegg.jp/pathway/{prefix}{pid}'
        urls.append(url)
    dic['url'] = urls
    dic.index.name = 'id'

    # Turns out names are unique too
    dic.to_csv(data_fdn + 'kegg_pathway_dataframe.tsv', sep='\t', index=True)
