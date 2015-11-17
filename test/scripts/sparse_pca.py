import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import decomposition, preprocessing
from itertools import combinations
from matplotlib import gridspec


# load data
exprs = np.log2(1 + pd.read_table(snakemake.input.est, index_col=[0, 1, 2])["expr_ev"].unstack(1).fillna(0))

# calculate PCA
spca = decomposition.PCA(n_components=3)
X = preprocessing.scale(exprs.transpose())
scores = pd.DataFrame(spca.fit_transform(X))
scores.columns = [
    "PC{} ({:.2%})".format(i + 1, expl_var)
    for i, expl_var in enumerate(spca.explained_variance_ratio_)
]
print(scores)
print(exprs)
scores.index = exprs.columns
pcs = list(combinations(scores.columns, 2))

# plot
sns.set_style("ticks")
fig = plt.figure(figsize=(9, 3))
gs = gridspec.GridSpec(1, 3)
for i, (a, b) in enumerate(pcs):
    ax = plt.subplot(gs[0, i])
    ax.plot(scores.loc[:, a], scores.loc[:, b], ".")
    #for group, df in annotation.groupby(annotation_column):
    #    ax.plot(scores.loc[df.index.values, a],
    #            scores.loc[df.index.values, b], ".",
    #            label=group)
    if i == 2:
        ax.legend(bbox_to_anchor=(1.6, 1))
    plt.xlabel(a)
    plt.ylabel(b)

# save the figure
sns.despine()
plt.tight_layout()
fig.savefig(snakemake.output[0], bbox_inches="tight")
