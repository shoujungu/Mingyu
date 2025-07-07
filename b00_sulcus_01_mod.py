import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
import pandas as pd
import numpy as np
import json
import scanpy as sc
import anndata as ad
import tensorflow.keras as keras

#-------------------------------------------------
n = 2000

fd_out = './out/b00_sulcus_01_mod'
f_in = f'./out/b00_sulcus_00_clus/data.h5ad'

##################################################
Path(fd_out).mkdir(parents=True, exist_ok=True)
ada = sc.read(f_in)

#-------------------------------------------------
class seq_dense:
    def __init__(self, n_layer=1, r=0.1, lr=0.0001, activation='swish', l1=0, l2=0):
        self.layer = n_layer
        #self.node=n_gene*2
        self.r = r
        self.lr = lr
        self.acti = activation
        self.l1 = l1
        self.l2 = l2
        return

    def build(self, n_in, n_out):
        #model
        reg = keras.regularizers.L1L2(l1=self.l1, l2=self.l2)
        model = keras.models.Sequential()
        model.add(keras.layers.Input((n_in,)))
        for i in range(self.layer):
            model.add(keras.layers.Dropout(self.r))
            model.add(keras.layers.Dense(n_in*2, activation=self.acti, kernel_initializer='he_normal', kernel_regularizer=reg))
        model.add(keras.layers.Dense(n_out, activation='softmax'))
        print(model.summary())
        #compile
        opt = keras.optimizers.Adam(learning_rate=self.lr)
        model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['categorical_accuracy'])
        return model

#-------------------------------------------------
#input genes, sort by expression ratio
df = ada.to_df().T
df['r'] = (df>0).mean(axis=1)
df = df.sort_values('r', ascending=False)
l_gene = df.index.tolist()[:n]

#training data
df = ada.to_df().reindex(l_gene, axis=1)
X = df.values

y = ada.obs['leiden'].astype('int')
y = keras.utils.to_categorical(y)

#model
n_out = len(ada.obs['leiden'].unique())
model = seq_dense()
model = model.build(n, n_out)
model.fit(X, y, epochs=5, batch_size=32, validation_split=0.1, callbacks=[])

#pred
l_col = [f'clus_{i}' for i in range(n_out)]
df_pred = pd.DataFrame(model.predict(X), columns=l_col, index=ada.obs.index)
df_pred.to_csv(f'{fd_out}/pred.csv')

print(df_pred)

