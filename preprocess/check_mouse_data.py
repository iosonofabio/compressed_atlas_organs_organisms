# vim: fdm=indent
'''
author:     Fabio Zanini
date:       09/12/22
content:    Check internal consistency of Emily's data before delving into the
            deep end.
'''
import sys
import pathlib
import h5py
import numpy as np
import pandas as pd
import anndata

import matplotlib.pyplot as plt


if __name__ == '__main__':

    if False:
        print('Check Emily\'s data')
        fdn_data = pathlib.Path('../data/Emilyslab')
        fn_atac = fdn_data / 'Mouse_Heart_ATAC_peaks.h5ad'
        fa = h5py.File(fn_atac)

        fn_rna = fdn_data / 'Mouse_Heart_RNA.h5ad'
        fr = h5py.File(fn_rna)

        age_a = fa['obs']['age'].asstr()[:]
        # FIXME: Y (young) -> 43 weeks -> 10 months, MA (middle age) -> 82 w -> 19m
        age_r = fr['obs']['age'].asstr()[:]

        cond_a = fa['obs']['condition'].asstr()[:]
        cond_r = fr['obs']['condition'].asstr()[:]

        # FIXME: ATAC "condition" is broken, just a copy of age
        # however, group combines age and condition, and is correct

        # RNA "condition" appears to be correct

        #with h5py.File(fn_atac) as fa, h5py.File(fn_rna) as fr:

        # Check cell annotations
        annos_a = fa['obs']['cellType'].asstr()[:]
        annos_r = fr['obs']['cellType'].asstr()[:]

        # FIX: RNA-Seq uses "Pericytes" instead of "Pericyte"
        annos_r[annos_r == 'Pericytes'] = "Pericyte"
        # FIX: RNA-Seq has a typo in lymphatics
        annos_r[annos_r == 'Lympathic_EC'] = "Lymphatic_EC"

        # FIX: ignore "unknowns"
        idx_a = ~pd.Index(annos_a).str.startswith("unknown")
        idx_r = ~pd.Index(annos_r).str.startswith("unknown")
        annos_a = annos_a[idx_a]
        annos_r = annos_r[idx_r]

        annos_a_count = pd.Series(annos_a).value_counts()
        annos_r_count = pd.Series(annos_r).value_counts()

        # Some renaming necessary
        # Neuronal -> Neuron
        # Coronary_EC -> Coronary
        # Lymphatic_EC -> Lymphatic
        # Myeloid -> Macrophage
        # Bcell -> B cell
        # Tcell -> T cell
        # Smooth_muscle -> Smooth muscle

        df_counts = pd.concat([annos_a_count, annos_r_count], axis=1).fillna(0).astype(int)
        df_counts.columns = ['ATAC', 'RNA']

        df_percent = 100 * df_counts / df_counts.sum()

        if False:
            fig, ax = plt.subplots()
            ax.scatter(0.01 + df_percent['ATAC'], 0.01 + df_percent['RNA'])
            ax.set_xlabel('ATAC [%]')
            ax.set_ylabel('RNA [%]')
            ax.set_xscale('log')
            ax.set_yscale('log')
            fig.tight_layout()
            plt.ion(); plt.show()

        # OK cell type percentages look ok

        # Now let's look at the matrices
        from scipy.sparse import csr_matrix

        if False:
            # ATAC matrix
            Xh5 = fa['raw']['X']
            X = csr_matrix(
                (np.asarray(Xh5['data']), np.asarray(Xh5['indices']), np.asarray(Xh5['indptr'])),
                )

            # FIXME: it's not binary, binarize ignoring diploidity
            print(X.max())
            X.data[:] = 1

            # FIXME: somehow it's float64, but we only need float32 for the averages
            # and really only bool (np.uint8) for the single cell data ;-)

        if False:
            # RNA matrix
            Xh5 = fr['raw']['X']
            X = csr_matrix(
                (np.asarray(Xh5['data']), np.asarray(Xh5['indices']), np.asarray(Xh5['indptr'])),
                )

            # FIXME: even though it says "raw", it's cptt -> natural logpp. Proof:
            X.data = np.exp(X.data) - 1
            print(X.sum(axis=1))
            nRNA = np.asarray(fr['obs']['nCount_RNA'])
            print(np.sort(((np.exp(X[0].data) - 1) * 1e-4) * nRNA[0])[-10:])

        if False:
            print('Check what "myeloid" cells actually are')
            Xh5 = fr['raw']['X']
            X = csr_matrix(
                (np.asarray(Xh5['data']), np.asarray(Xh5['indices']), np.asarray(Xh5['indptr'])),
                )
            genes = fr['raw']['var']['_index'].asstr()[:]
            cells = fr['obs']['_index'].asstr()[:]
            celltypes = fr['obs']['cellType'].asstr()[:]
            adata = anndata.AnnData(X=X)
            adata.obs_names = cells
            # The last 3 genes are never seen, so the sparse matrix constructor
            # infers it wrong. No biggie... I think
            adata.var_names = genes[:X.shape[1]]
            adata.obs['cellType'] = celltypes

            adata_myeloid = adata[adata.obs['cellType'] == 'Myeloid']
            markers = ['Ptprc', 'Cd14', 'Marco', 'Mrc1', 'Ms4a7',
                       'C1qa', 'C1qb', 'C1qc', 'Ccl12', 'Ccl24']
            print((adata_myeloid[:, ].X > 0).mean(axis=0))

            # They are interstitial macrophages

        if True:
            print('Compute markers for vent VS atrial cardiomyocytes')
            Xh5 = fr['raw']['X']
            X = csr_matrix(
                (np.asarray(Xh5['data']), np.asarray(Xh5['indices']), np.asarray(Xh5['indptr'])),
                )
            genes = fr['raw']['var']['_index'].asstr()[:]
            cells = fr['obs']['_index'].asstr()[:]
            celltypes = fr['obs']['cellType'].asstr()[:]
            adata = anndata.AnnData(X=X)
            adata.obs_names = cells
            # The last 3 genes are never seen, so the sparse matrix constructor
            # infers it wrong. No biggie... I think
            adata.var_names = genes[:X.shape[1]]
            adata.obs['cellType'] = celltypes

            adata_vent = adata[adata.obs['cellType'] == 'Vent']
            adata_atr = adata[adata.obs['cellType'] == 'Atrial']

            avg_vent = pd.Series(
                    np.asarray(adata_vent.X.mean(axis=0))[0],
                    index=adata.var_names)
            avg_atr = pd.Series(
                    np.asarray(adata_atr.X.mean(axis=0))[0],
                    index=adata.var_names)
            diff = (avg_atr - avg_vent)
            #https://link.springer.com/protocol/10.1007/978-1-0716-1484-6_14
            markers_vent = ['Hey2', 'Irx4', 'Myl2']
            markers_atr = ['Nr2f1', 'Nr2f2', 'Nppa', 'Tbx5', 'Kcna5', 'Kcnj3']
            print(diff.nlargest(10)) # Lrp1b, Gsn, Ebf1, Dcn, Sec61a1, Pcdh9, Apoe
            print(diff.nsmallest(10)) # Lrrtm3, Fgfr2, Mhrt, Trdn, Myh7

    if True:
        print('Check TMS data consistency with Emily\'s')
        fdn_data = pathlib.Path('../data/tabula_muris_senis')
        fn_rna = fdn_data / 'tabula-muris-senis-droplet-processed-official-annotations-Heart_and_Aorta.h5ad'
        adata = anndata.read(fn_rna)
