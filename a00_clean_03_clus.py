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
res = 0.5
sample = '3M'

fd_out = './out/a00_clean_03_clus'
f_in = './out/a00_clean_02_filter/data.h5ad'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
ada = sc.read(f_in)

#-------------------------------------------------
def plt_umap(ada, f_out, title=None, cmap='tab20', sz=(6, 5), s=5, obsm='X_umap', col='leiden', ncols=1, l_box=[1, 0.8], leg_title='Clusters', leg=True):
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
#clus
sc.pp.neighbors(ada)
sc.tl.leiden(ada, resolution=res)
ada.write(f'{fd_out}/data.h5ad')

#plot
f_out = f'{fd_out}/clus.png'
plt_umap(ada, f_out, title=sample)


#plot clus
Path(f'{fd_out}/clus').mkdir(parents=True, exist_ok=True)
ada.obs['leiden'] = ada.obs['leiden'].astype('str')
for i in ada.obs['leiden'].unique():
    #prep
    ada.obs['tmp'] = 'Others'
    ada.obs.loc[ada.obs['leiden']==i, ['tmp']] = f'Cluster {i}'
    ada.obs['tmp'] = pd.Categorical(ada.obs['tmp'], categories=[f'Cluster {i}', 'Others'], ordered=True) 
    #plot
    f_out = f'{fd_out}/clus/clus_{i}.png'
    title = f'Cluster {i}'
    plt_umap(ada, f_out, title=title, cmap=['#ff0000', '#888888'], leg=False, col='tmp')
