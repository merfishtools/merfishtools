import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


sns.set(style="ticks", context=snakemake.wildcards.context)
sns.set_palette("Dark2", n_colors=7)

x, y = snakemake.config["plots"]["figsize"]
plt.figure(figsize=(x, y * 2))

exprs = [pd.read_table(f, index_col=0).stack(0) for f in snakemake.input.exprs]
exprs = pd.concat(exprs, keys=snakemake.params.expmnts)

diffexp = pd.read_table(snakemake.input.diffexp, index_col=0)
diffexp.sort_values("cv_ev", inplace=True)

significant = diffexp[diffexp["diff_fdr"] <= 0.05]
significant = significant[~significant.index.str.startswith("blank") & ~significant.index.str.startswith("notarget")]

exprs = exprs.loc[:, significant.index]

def subplot(exprs, xlim):
    means = exprs.groupby(level=[0, 1]).mean()

    #high_expmnts = means.groupby(level=[0, 1]).max().index
    #exprs["col"] = "black"
    #exprs.loc[high_expmnts, "col"] = "red"
    exprs = exprs.reset_index()
    exprs.columns = ["expmnt", "feat", "cell", "expr_ev"]#, "col"]
    exprs["col"] = "grey"
    #exprs["col"] = exprs["expmnt"]
    #exprs.loc[exprs["expr_ev"] < 10, "col"] = np.nan
    means = means.reset_index()
    means.columns = ["expmnt", "feat", "mean"]

    #exprs.reindex(np.random.permutation(exprs.index))

    #plt.subplot2grid((6, 1), (i, 0), rowspan=5 if i > 0 else 1)
    ax = sns.stripplot(x="expr_ev", y="feat", hue="col", data=exprs, alpha=0.6, jitter=True, size=1, clip_on=True, palette=sns.color_palette("Greys", n_colors=1))
    ax = sns.stripplot(x="mean", y="feat", hue="expmnt", data=means, ax=ax, zorder=5, size=5, jitter=True, clip_on=True)#, palette=sns.color_palette("Set1", desat=0.5))
    #ax.set_xscale("log")
    ax.legend().set_visible(False)
    plt.xlabel("mean expression")
    plt.ylabel("")
    sns.despine()
    plt.xlim(xlim)

means = exprs.groupby(level=1).mean()
high = means[means > 200].index
low = means[means <= 200].index

#subplot(exprs.loc[:, high], 1, (200, 2000))
subplot(exprs.loc[:, low], (0, 50))
plt.savefig(snakemake.output[0], bbox_inches="tight")
