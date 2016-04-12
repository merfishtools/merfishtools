import pandas as pd
import matplotlib.pyplot as plt

from merfishtools.estimates import plot_estimates


def read_cdf(path_or_buffer):
    cdf = pd.read_table(path_or_buffer)
    # set all but the last two values as index
    cdf.set_index(cdf.columns[:-2], inplace=True)
    return cdf


def plot_cdf(cdf, expected_value=None, credible_interval=None, legend=True):
    plt.step(cdf.iloc[:, 0], cdf.iloc[:, 1], "k-", label="", clip_on=False, zorder=6)
    ylim = plt.ylim()
    plot_estimate(ylim, expected_value, credible_interval, legend=legend)
    if credible_interval is not None or expected_value is not None:
        mask_estimate(ylim, cdf)
    plt.ylabel("CDF")


def plot_pmf(cdf, expected_value=None, credible_interval=None, legend=True):
    pmf = cdf
    pmf.iloc[:, 1][1:] = pmf.iloc[:, 1][1:] - pmf.iloc[:, 1][:-1]
    plt.plot(pmf.iloc[:, 0], pmf.iloc[:, 1], "ko", label="", ms=4, clip_on=False, zorder=6)
    plot_estimate(plt.ylim(), expected_value, credible_interval, legend=legend)
    if credible_interval is not None or expected_value is not None:
        mask_estimate(ylim, cdf)
    plt.ylabel("PMF")


def mask_estimate(ylim, cdf):
    plt.fill_between(cdf.iloc[: 0], cdf.iloc[: 1], ylim[1],  zorder=5, facecolor="white", edgecolor="white", step="pre")
    plt.ylim(ylim)
