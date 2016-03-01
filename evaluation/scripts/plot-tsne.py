import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing, manifold
from itertools import combinations
from matplotlib import gridspec


# load data
exprs = [pd.read_table(f, index_col=0) for f in snakemake.input]
exprs = pd.concat(exprs, axis="columns", keys=range(1, len(snakemake.input) + 1))
exprs = np.log10(1 + exprs.transpose())

# calculate t-SNE embedding
tsne = manifold.TSNE(random_state=21498)
X = preprocessing.scale(exprs)
embedding = pd.DataFrame(tsne.fit_transform(X), columns=["x", "y"])
embedding.index = exprs.index

# plot
sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
width, height = snakemake.config["plots"]["figsize"]
fig = plt.figure(figsize=snakemake.config["plots"]["figsize"])

markers = dict(zip(snakemake.params.codebooks, "o^"))
for color, codebook, (expmnt, _embedding) in zip(sns.color_palette("muted", len(snakemake.params.codebooks)),
                                              snakemake.params.codebooks, embedding.groupby(level=0)):
    plt.scatter(_embedding['x'], _embedding['y'],
               marker=markers[codebook], label=expmnt,
               c=color, edgecolors="face", alpha=0.7)


# save the figure
plt.axis("off")
extent = plt.gca().get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
#plt.tight_layout()
plt.savefig(snakemake.output[0], bbox_inches=extent, pad_inches=0)
