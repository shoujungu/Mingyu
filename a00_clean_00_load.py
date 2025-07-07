import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
import pandas as pd
import numpy as np
import json
import scanpy as sc
import anndata as ad
from multiprocessing import Pool, Manager

#-------------------------------------------------
sample = '3M'
idx = 'GSM8446833'

fd_out = './out/a00_clean_00_load'
fd_in = f'../a00_raw/data/{sample}'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)

#-------------------------------------------------
#load
prefix = f'{idx}_{sample}_Cochlea_'
ada = sc.read_10x_mtx(fd_in, cache=False, prefix=prefix)
sc.pp.filter_genes(ada, min_counts=1)

#metrics
ada.var['mt'] = ada.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(ada, qc_vars=['mt'], percent_top=None, inplace=True, log1p=False)

#clean obs
ada.obs['sample'] = sample
ada.obs.index = ada.obs.index.map(lambda x: f'{x}_{sample}')
ada.obs = ada.obs.rename(columns={'n_genes_by_counts': 'n_gene', 'total_counts': 'cnt', 'total_counts_mt': 'n_mt', 'pct_counts_mt': 'pct_mt'})
ada.obs['pct_mt'] = ada.obs['pct_mt'].round(2)
ada.obs['cnt'] = ada.obs['cnt'].astype('int')
ada.obs['n_mt'] = ada.obs['n_mt'].astype('int')

#clean var
ada.var = ada.var.drop(['feature_types', 'mean_counts', 'total_counts'], axis=1)
ada.var = ada.var.rename(columns={'gene_ids': 'ids', 'n_counts': 'cnt', 'n_cells_by_counts': 'n_cell', 'pct_dropout_by_counts': 'pct_drop'})
ada.var['cnt'] = ada.var['cnt'].astype('int')
ada.var['pct_drop'] = ada.var['pct_drop'].round(2)

#save
print(ada.obs)
ada.write(f'{fd_out}/data.h5ad') 

