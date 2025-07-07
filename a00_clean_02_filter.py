import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
import pandas as pd
import numpy as np
import json
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool, Manager
from scipy.stats import pearsonr, zscore

#-------------------------------------------------
fd_out = './out/a00_clean_02_filter'
f_in = './out/a00_clean_01_qc/data.h5ad'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
ada = sc.read(f_in)
print(ada)

#-------------------------------------------------
#norm
sc.pp.normalize_total(ada, exclude_highly_expressed=True)
sc.pp.log1p(ada)

#filter genes
ada = ada[:, ada.var['filter']==0]

#filter cells
ada = ada[ada.obs['filter_cnt']==0, :]
ada = ada[ada.obs['filter_gene']==0, :]
ada = ada[ada.obs['filter_mt']==0, :]
ada = ada[ada.obs['filter_dd']==0, :].copy()

#clean
ada.var = ada.var.reindex(['ids'], axis=1)
ada.obs = ada.obs.reindex(['sample'], axis=1)

#PCA
ada_c = ada.copy()
sc.pp.highly_variable_genes(ada_c, n_top_genes=2000, inplace=True)
ada_c = ada_c[:, ada_c.var.highly_variable]
sc.tl.pca(ada_c, svd_solver='arpack', n_comps=30)
ada.obsm['X_pca'] = ada_c.obsm['X_pca'] 

#UMAP
sc.pp.neighbors(ada_c)
sc.tl.umap(ada_c, n_components=2, random_state=42) 
ada.obsm['X_umap'] = ada_c.obsm['X_umap']

#write
print(ada)
ada.write(f'{fd_out}/data.h5ad')
