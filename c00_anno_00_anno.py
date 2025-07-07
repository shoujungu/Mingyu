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

#-------------------------------------------------
sample = '3M'
l_sub = [
'./out/b00_sulcus_04_anno/data.h5ad',
'./out/b01_unknown_01_anno/data.h5ad'
]

fd_out = './out/c00_anno_00_anno'
f_ada = './out/a00_clean_09_anno/data.h5ad'
f_meta = '../a00_raw/meta/cell_3m.json'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
with open(f_meta, 'r') as f: d_meta = json.load(f)
ada = sc.read(f_ada)

#-------------------------------------------------
def plt_umap(ada, f_out, title=None, cmap='tab20', sz=(6, 5), s=5, obsm='X_umap', col='cell', ncols=1, l_box=[1, 0.8], leg_title='', leg=True):
    #prep
    df = pd.DataFrame(ada.obsm[obsm], columns=['x', 'y'], index=ada.obs.index)
    df[col] = ada.obs[col]
    #plot
    fig, ax = plt.subplots(figsize=sz)
    ax=sns.scatterplot(x='x', y='y', data=df, hue=col, palette=cmap, s=s, linewidth=0, alpha=0.7)
    #adjust
    plt.axis('off')
    ax.set_title(title, fontsize=24, weight='semibold', pad=12)
    #legend
#    reorder = lambda l, nc: sum((l[i::nc] for i in range(nc)), [])
#    h, l = ax.get_legend_handles_labels()
    if leg:
        plt.legend(loc='upper left', bbox_to_anchor=l_box, ncols=ncols, title=leg_title, frameon=True, prop={'size': 8, 'weight': 'semibold'}, markerscale=3, title_fontproperties={'size': 10, 'weight': 'semibold'})
    else:
        plt.legend().remove()
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#-------------------------------------------------
#anno
ada.obs['cell'] = ada.obs['cell'].astype('str')
for fname in l_sub:
    ada_sub = sc.read(fname)
    ada_sub.obs['cell'] = ada_sub.obs['cell'].astype('str')
    ada.obs.loc[ada_sub.obs.index, ['cell']] = ada_sub.obs['cell']

ada.obs['cell'] = pd.Categorical(ada.obs['cell'], categories=list(d_meta), ordered=True)
ada.write(f'{fd_out}/data.h5ad')

#plot
f_out = f'{fd_out}/cell.png'
cmap = list(d_meta.values())
plt_umap(ada, f_out, title=sample, cmap=cmap)

#-------------------------------------------------
ada = sc.read(f'{fd_out}/data.h5ad')
print(ada.obs['cell'].astype('str').unique().tolist())

#plot clus
Path(f'{fd_out}/cell').mkdir(parents=True, exist_ok=True)
ada.obs['cell'] = ada.obs['cell'].astype('str')
for cell in ada.obs['cell'].unique():
    #prep
    ada.obs['tmp'] = 'Others'
    ada.obs.loc[ada.obs['cell']==cell, ['tmp']] = cell
    ada.obs['tmp'] = pd.Categorical(ada.obs['tmp'], categories=[cell, 'Others'], ordered=True) 
    #plot
    f_out = f'{fd_out}/cell/{cell}.png'
    title = cell
    plt_umap(ada, f_out, title=title, cmap=['#ff0000', '#888888'], leg=False, col='tmp')
