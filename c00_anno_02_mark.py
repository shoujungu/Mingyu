from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
from multiprocessing import Pool

#-------------------------------------------------
fd_out = './out/c00_anno_02_mark'
f_ada = './out/c00_anno_00_anno/data.h5ad'
f_pred = './out/c00_anno_01_mod/pred.csv'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)

#load
ada = sc.read(f_ada)
df = ada.to_df()
df_pred = pd.read_csv(f_pred, index_col=0)

#-------------------------------------------------
def mainf(cell):
    print(cell)
    #corr
    df_corr = pd.DataFrame(df.corrwith(df_pred[cell], axis=0), columns=['corr'])
    df_corr = df_corr.sort_values('corr', ascending=False)
    df_corr.to_csv(f'{fd_out}/{cell}.csv')
    return

#-------------------------------------------------
l_cell = ada.obs['cell'].astype('str').unique().tolist()
with Pool(8) as p: p.map(mainf, l_cell)

