#%%
import glob
import pandas as pd
import scanpy as sc
import os
#%%
def read_file(GSE):
    file = glob.glob(f'Data/{GSE}/*_raw_UMI_count_matrix.txt.gz')[0]
    adata = sc.read_text(file, delimiter='\t').T
    print(f'Loaded {adata.shape[0]} cells and {adata.shape[1]} genes')
    return adata

def normalize(adata, target_sum=10000, n_top_genes=1000):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    # sc.pl.highly_variable_genes(adata)

def pca(adata, n_comps=None):
    sc.tl.pca(adata, n_comps=n_comps)
    # sc.pl.pca_variance_ratio(adata, n_pcs=min(50, n_comps))

def magic(GSE):
    adata = read_file(GSE)
    normalize(adata)
    pca(adata)
    sc.external.pp.magic(adata, n_jobs=os.cpu_count()-2)
    # Save the results
    # adata.write(f'Data/{GSE}-magic.h5ad')
    os.mkdir(f'Data/{GSE}-magic')
    adata.write_csvs(f'Data/{GSE}-magic/', skip_data=False)

GSEs = ['GSE132257']
for GSE in GSEs:
    magic(GSE)