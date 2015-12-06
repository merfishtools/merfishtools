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
exprs = [pd.read_table(snakemake.input[0], index_col=0) for f in snakemake.input]
exprs = pd.concat(exprs, axis="columns", keys=range(1, len(snakemake.input) + 1))
exprs = np.log10(1 + exprs.transpose())

# calculate PCA
spca = decomposition.PCA(n_components=3)
X = preprocessing.scale(exprs)
scores = pd.DataFrame(spca.fit_transform(X))
scores.columns = [
    "PC{} ({:.2%})".format(i + 1, expl_var)
    for i, expl_var in enumerate(spca.explained_variance_ratio_)
]
scores.index = exprs.index
pcs = list(combinations(scores.columns, 2))

# plot
sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
width, height = snakemake.config["plots"]["figsize"]
plt.figure(figsize=(width * 3, height))

gs = gridspec.GridSpec(1, 3)
for i, (a, b) in enumerate(pcs):
    ax = plt.subplot(gs[0, i])
    for expmnt, _scores in scores.groupby(level=0):
        ax.plot(_scores.loc[:, a], _scores.loc[:, b], ".", label=expmnt, markersize=4, alpha=0.9)
    #if i == 2:
    #    ax.legend(bbox_to_anchor=(1.6, 1), title="experiment")
    plt.xlabel(a)
    plt.ylabel(b)

# save the figure
sns.despine()
plt.tight_layout()
plt.savefig(snakemake.output[0], bbox_inches="tight")
