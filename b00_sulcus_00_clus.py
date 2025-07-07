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

#-------------------------------------------------
cell = 'Sulcus'
res = 0.2

fd_out = './out/b00_sulcus_00_clus'
f_in = './out/a00_clean_09_anno/data.h5ad'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
ada_all = sc.read(f_in)

#-------------------------------------------------
def plt_umap(ada, f_out, title=None, cmap='tab20', sz=(6, 5), s=20, obsm='X_umap', col='leiden', ncols=1, l_box=[1, 0.8], leg_title='Clusters', leg=True):
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
#prep
ada = ada_all[ada_all.obs['cell']==cell, :].copy()
ada.obs = ada.obs.drop('cell', axis=1)

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

#clus
sc.pp.neighbors(ada)
sc.tl.leiden(ada, resolution=res)
ada.write(f'{fd_out}/data.h5ad')

#plot
f_out = f'{fd_out}/sub_clus.png'
plt_umap(ada, f_out, title='Sub-Cluster', leg=False)

#-------------------------------------------------
#all clus
for clus in ada.obs['leiden'].unique().tolist():
    ada_all.obs['tmp'] = 'Others'
    df_clus = ada.obs.query('leiden==@clus')
    ada_all.obs.loc[df_clus.index, ['tmp']] = 'Clus'
    ada_all.obs['tmp'] = pd.Categorical(ada_all.obs['tmp'], categories=['Clus', 'Others'], ordered=True)   
    #plot
    f_out = f'{fd_out}/{clus}.png'
    plt_umap(ada_all, f_out, title=clus, leg=False, s=5, cmap=['#ff0000', '#888888'], col='tmp')


