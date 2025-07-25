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
import matplotlib.colors as mcolors

#-------------------------------------------------
n = 5

fd_out = './out/c00_anno_04_plt-exp'
f_in = f'./out/c00_anno_00_anno/data.h5ad'
fd_gene = './out/c00_anno_02_mark'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
l_fname = list(Path(fd_gene).glob('*.csv'))
ada = sc.read(f_in)

#prep
df_xy = pd.DataFrame(ada.obsm['X_umap'], columns=['x', 'y'], index=ada.obs.index)
df = ada.to_df()
df = df.merge(df_xy, left_index=True, right_index=True)

#-------------------------------------------------
def truncate_cmap(cmap, min_v=0.25, max_v=1.0, n=100):
    #from ChatGPT
    cmap = plt.get_cmap(cmap)
    new_cmap = mcolors.LinearSegmentedColormap.from_list('name', cmap(np.linspace(min_v, max_v, n)))
    return new_cmap

def plt_exp(gene, f_out, df=df, sz=(3, 2.8), s=3, cmap='BuPu', vmax=None):
    cmap = truncate_cmap(cmap)
    #prep 
    if vmax is None: vmax = max(int(df[gene].max()), 3)
    norm = plt.Normalize(0, vmax)
    l_tick=[0, vmax]
    #plot  (left, bottom, width, height)
    fig = plt.figure(figsize=sz)
    ax = fig.add_axes([0, 0, 0.85, 0.9])
    ax = sns.scatterplot(x='x', y='y', data=df, hue=gene, palette=cmap, s=s, linewidth=0, alpha=0.9, hue_norm=norm)
    #adjust
    plt.axis('off')
    ax.set_title(gene, fontsize=16, weight='semibold', style='italic')
    ax.get_legend().remove()
    #color bar (https://stackoverflow.com/questions/62884183/trying-to-add-a-colorbar-to-a-seaborn-scatterplot)
    ax_cbar = fig.add_axes([0.85, 0.3, 0.02, 0.4])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, cax=ax_cbar, ticks=l_tick)
    cb.ax.tick_params(labelsize=11)
    #save
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

def mainf(fname):
    clus = Path(fname).stem
    print(clus)
    df = pd.read_csv(fname, index_col=0)
    l_gene = df.index.tolist()[:n]
    #plot
    Path(f'{fd_out}/{clus}').mkdir(parents=True, exist_ok=True)
    for gene in l_gene:
        f_out = f'{fd_out}/{clus}/{gene}.png'
        plt_exp(gene, f_out)
    return

#-------------------------------------------------
with Pool(10) as p: p.map(mainf, l_fname)
