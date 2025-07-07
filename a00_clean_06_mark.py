from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
from multiprocessing import Pool

#-------------------------------------------------
fd_out = './out/a00_clean_06_mark'
f_ada = './out/a00_clean_03_clus/data.h5ad'
f_pred = './out/a00_clean_05_mod/pred.csv'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)

#load
ada = sc.read(f_ada)
df = ada.to_df()
df_pred = pd.read_csv(f_pred, index_col=0)

#-------------------------------------------------
def mainf(clus):
    print(clus)
    #corr
    df_corr = pd.DataFrame(df.corrwith(df_pred[f'clus_{clus}'], axis=0), columns=['corr'])
    df_corr = df_corr.sort_values('corr', ascending=False)
    df_corr.to_csv(f'{fd_out}/clus_{clus}.csv')
    return

#-------------------------------------------------
l_clus = ada.obs['leiden'].astype('str').unique().tolist()
with Pool(8) as p: p.map(mainf, l_clus)

