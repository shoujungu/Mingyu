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
'0': 'Sulcus-1',
'1': 'Sulcus-2',
'2': 'Sulcus-3'
}

fd_out = './out/b00_sulcus_04_anno'
f_in = './out/b00_sulcus_00_clus/data.h5ad'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
ada = sc.read(f_in)

#-------------------------------------------------
#anno
ada.obs['cell'] = ada.obs['leiden'].apply(lambda x: d_anno[x])
ada.obs = ada.obs.drop('leiden', axis=1)

print(ada.obs)
ada.write(f'{fd_out}/data.h5ad')
