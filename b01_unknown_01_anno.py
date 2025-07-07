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
import json

#-------------------------------------------------
d_anno = {
'0': 'Unknown',
'1': 'Fibrocyte',
'2': 'SVB'
}

fd_out = './out/b01_unknown_01_anno'
f_in = './out/b01_unknown_00_clus/data.h5ad'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
ada = sc.read(f_in)

#-------------------------------------------------
#anno
ada.obs['cell'] = ada.obs['leiden'].apply(lambda x: d_anno[x])
ada.obs = ada.obs.drop('leiden', axis=1)

print(ada.obs)
ada.write(f'{fd_out}/data.h5ad')
