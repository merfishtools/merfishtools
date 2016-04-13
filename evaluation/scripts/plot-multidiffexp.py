import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid.inset_locator import inset_axes


sns.set(style="ticks", context=snakemake.wildcards.context)
sns.set_palette("Dark2", n_colors=7)

fig_x, fig_y = snakemake.config["plots"]["figsize"]
figsize = (fig_x * 4, fig_y * 4)
plt.figure(figsize=figsize)

exprs = [pd.read_table(f, index_col=0).stack(0) for f in snakemake.input.exprs]
exprs = pd.concat(exprs, keys=snakemake.params.expmnts)

diffexp = pd.read_table(snakemake.input.diffexp, index_col=0)
diffexp.sort_values("cv_ev", inplace=True)

pmfs = pd.read_table(snakemake.input.diffexp_pmf, index_col=0)

significant = diffexp[diffexp["diff_fdr"] <= 0.05]
significant = significant[~significant.index.str.startswith("blank") & ~significant.index.str.startswith("notarget")]

exprs = exprs.loc[:, significant.index]

n = exprs.index.levels[1].size
for i, (gene, gene_exprs) in enumerate(exprs.groupby(level=1)):
    vmax = gene_exprs.quantile(0.95)
    ax = plt.subplot(4, 4, i + 1)
    for _, exp_exprs in gene_exprs.groupby(level=0):
        sns.kdeplot(exp_exprs, ax=ax, linewidth=1, alpha=0.5)
    if i == n - 1:
        plt.xlabel("gene expression")
    plt.ylabel(gene)
    plt.xlim((0, vmax))
    sns.despine()
    plt.yticks([])

    # plot cdf
    cdf_ax = inset_axes(ax, width="30%", height=0.5, loc=1)
    pmf = pmfs.loc[gene]
    cdf = np.cumsum(np.exp(pmf["prob"]))
    plt.step(pmf["cv"], cdf, "k-", clip_on=False, zorder=6)
    ylim = plt.ylim()
    est = significant.loc[gene]

    ev = est["cv_ev"]
    ci_lower, ci_upper = est[["cv_ci_lower", "cv_ci_upper"]]

    plt.fill([ci_lower, ci_upper, ci_upper, ci_lower], [0, 0, ylim[1], ylim[1]], "red", lw=0, label="95% credible interval", alpha=0.5)
    plt.vlines([ev], *ylim, colors="red", linestyles="-", label="expected value")

    plt.fill_between(pmf["cv"], cdf, 1.2,  zorder=5, facecolor="white", edgecolor="white", step="pre")
    plt.ylim(ylim)  
    #plt.xlim((0, plt.xlim()[1]))
    plt.setp(cdf_ax.get_xticklabels(), rotation=45, ha="right")
    cdf_ax.tick_params(pad=1)
    plt.locator_params(nbins=4)
    sns.despine()


plt.tight_layout()
plt.savefig(snakemake.output[0], bbox_inches="tight")

"""
plt.figure(figsize=figsize)
pmfs = pmfs.loc[significant.index]
for i, (gene, pmf) in enumerate(pmfs.groupby(level=0)):
    plt.subplot(4, 4, i + 1)
    cdf = np.cumsum(np.exp(pmf["prob"]))
    plt.step(pmf["cv"], cdf, "k-", clip_on=True, zorder=6)
    plt.xlabel(gene)
    ylim = plt.ylim()
    est = significant.loc[gene]

    ev = est["cv_ev"]
    ci_lower, ci_upper = est[["cv_ci_lower", "cv_ci_upper"]]

    plt.fill([ci_lower, ci_upper, ci_upper, ci_lower], [0, 0, ylim[1], ylim[1]], "red", lw=0, label="95% credible interval", alpha=0.5)
    plt.vlines([ev], *ylim, colors="red", linestyles="-", label="expected value")

    plt.fill_between(pmf["cv"], cdf, 1.2,  zorder=5, facecolor="white", edgecolor="white", step="pre")
    plt.ylim(ylim)
    #plt.xlim((0, plt.xlim()[1]))
    sns.despine()
plt.tight_layout()

plt.savefig(snakemake.output.cdf, bbox_inches="tight")"""
