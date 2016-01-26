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
exprs = [pd.read_table(f, index_col=0) for f in snakemake.input]
exprs = pd.concat(exprs, axis="columns", keys=range(1, len(snakemake.input) + 1))
exprs = np.log10(1 + exprs.transpose())

# calculate PCA
spca = decomposition.PCA(n_components=3)
X = preprocessing.scale(exprs)
scores = pd.DataFrame(spca.fit_transform(X))
# scores.columns = ["PC{}".format(i) for i in range(1, 4)]
scores.columns = [
    "PC{} ({:.2%})".format(i + 1, expl_var)
    for i, expl_var in enumerate(spca.explained_variance_ratio_)
]
scores.index = exprs.index
pcs = list(combinations(scores.columns, 2))

# plot
sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
width, height = snakemake.config["plots"]["figsize"]
plt.figure(figsize=(height * 3, height))

markers = dict(zip(snakemake.params.codebooks, "o^"))
print(markers)
gs = gridspec.GridSpec(1, 3)
for i, (a, b) in enumerate(pcs):
    ax = plt.subplot(gs[0, i])
    for color, codebook, (expmnt, _scores) in zip(sns.color_palette("muted", len(snakemake.params.codebooks)),
                                                  snakemake.params.codebooks, scores.groupby(level=0)):
        ax.scatter(_scores.loc[:, a], _scores.loc[:, b],
                   marker=markers[codebook], label=expmnt,
                   c=color, edgecolors="face", alpha=0.7)
    #if i == 2:
    #    ax.legend(bbox_to_anchor=(1.6, 1), title="experiment")
    plt.xlabel(a)
    plt.ylabel(b)

# save the figure
sns.despine()
plt.tight_layout()
plt.savefig(snakemake.output[0], bbox_inches="tight")
