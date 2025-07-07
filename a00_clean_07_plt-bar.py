import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool, Manager

#-------------------------------------------------
n = 20

fd_out = './out/a00_clean_07_plt-bar'
fd_in = './out/a00_clean_06_mark'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
l_fname = list(Path(fd_in).glob('*.csv'))

#-------------------------------------------------
def plt_bar(df, f_out, title=None, sz=(2.8, 5), clr=None):
    #plot
    fig, ax = plt.subplots(figsize=sz)
    ax = sns.barplot(data=df, x='corr', y=df.index, color='grey')
    
    #adjust
    plt.title(title, fontsize=22, weight='semibold', pad=12, x=0.45)
    plt.xlabel('Correlation', fontsize=12, weight='semibold')
    plt.ylabel('', fontsize=12, weight='semibold')
    plt.xlim([0, 0.8])
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8], fontsize=8)
    plt.yticks(fontsize=9.5, weight='semibold')
    
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#-------------------------------------------------
for fname in l_fname:
    clus = Path(fname).stem
    df = pd.read_csv(fname, index_col=0)
    df = df.iloc[:n, :]

    #plot
    f_out = f'{fd_out}/{clus}.png'
    plt_bar(df, f_out, title=f'Cluster {clus.split("_")[-1]}')
