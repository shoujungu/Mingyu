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
'0': 'TBC',
'1': 'TBC',
'2': 'Sulcus',
'3': 'Unknown',
'4': 'IHC',
'5': 'Interdental',
'6': 'Reissner',
'7': 'OHC',
'8': 'SGN',
'9': 'Unknown',
'10': 'Support',
'11': 'Macrophage',
'12': 'Schwann',
'13': 'Endothelial'
}

fd_out = './out/a00_clean_09_anno'
f_in = './out/a00_clean_03_clus/data.h5ad'
f_meta = '../a00_raw/meta/cell.json'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
with open(f_meta, 'r') as f: d_meta = json.load(f)
ada = sc.read(f_in)

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
ada.obs['cell'] = ada.obs['leiden'].apply(lambda x: d_anno[x])
ada.obs = ada.obs.drop('leiden', axis=1)

#cat
l_cell = [i for i in list(d_meta) if i in list(d_anno.values())]
ada.obs['cell'] = pd.Categorical(ada.obs['cell'], categories=l_cell, ordered=True)
ada.write(f'{fd_out}/data.h5ad')

#plot
f_out = f'{fd_out}/cell.png'
cmap = [d_meta[i] for i in l_cell]
plt_umap(ada, f_out, title='Cell', cmap=cmap)
