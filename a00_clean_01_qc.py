import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
import pandas as pd
import numpy as np
import json
import scanpy as sc
import anndata as ad
from multiprocessing import Pool, Manager
import re

#-------------------------------------------------
fd_out = './out/a00_clean_01_qc'
f_in = './out/a00_clean_00_load/data.h5ad'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
ada = sc.read(f_in)

#-------------------------------------------------
def filter_gene(l_gene, l_ignore=['Ercc-', 'Tg-', 'Tg_', 'n-', '-ps', '-as', '.']):
    #general
    for t in l_ignore: l_gene=[i for i in l_gene if not t in i]  
    l_gene=[i for i in l_gene if not re.match('^Gm[0-9][0-9].*', i)]
    l_gene=[i for i in l_gene if not re.match('^[a-z].*', i)]
    l_gene=[i for i in l_gene if not re.match('^[A-Z][A-Z].*', i)]
    l_gene=[i for i in l_gene if len(i)<15] 
    #other
    l_gene=[i for i in l_gene if not re.match('^Rp[l,s].*', i)]
    l_gene=[i for i in l_gene if not re.match('^Mir[0-9][0-9].*', i)]
    l_gene=[i for i in l_gene if not re.match('^Linc[0-9][0-9].*', i)]
    l_gene=[i for i in l_gene if not (('Rik' in i) & (not i.startswith('Rik')))]
    l_gene=[i for i in l_gene if not 'mt-' in i]
    #sort
    l_gene.sort()
    return l_gene

def mainf(ada, min_cnt=100, max_cnt=10000, min_gene=100, max_gene=5000, max_mt=10):
    #qc cell (0 is keep, 1 is filter)
    ada.obs['filter_cnt'] = 0
    ada.obs.loc[(ada.obs['cnt']<min_cnt) | (ada.obs['cnt']>max_cnt), ['filter_cnt']] = 1
    ada.obs['filter_gene'] = 0
    ada.obs.loc[(ada.obs['n_gene']<min_gene) | (ada.obs['n_gene']>max_gene), ['filter_gene']] = 1
    ada.obs['filter_mt'] = 0
    ada.obs.loc[ada.obs['pct_mt']>max_mt, ['filter_mt']] = 1
    #dedoublet
    sc.pp.scrublet(ada)
    ada.obs['filter_dd'] = ada.obs['predicted_doublet'].astype('int')
    ada.obs = ada.obs.drop('predicted_doublet', axis=1)
    #qc gene (0 is keep, 1 is filter)
    ada.var['filter'] = 0
    l_keep = filter_gene(ada.var.index.tolist())
    ada.var.loc[~ada.var.index.isin(l_keep), ['filter']] = 1
    #save
    ada.write(f'{fd_out}/data.h5ad')
    return

#-------------------------------------------------
mainf(ada)
