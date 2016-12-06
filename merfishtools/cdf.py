import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from merfishtools.estimates import plot_estimate


def read_cdf(path_or_buffer):
    cdf = pd.read_table(path_or_buffer)
    # set all but the last two values as index
    cdf.set_index(list(cdf.columns[:-2]), inplace=True)
    cdf.sort_index(inplace=True, kind="mergesort")
    return cdf


def plot_cdf(cdf, expected_value=None, credible_interval=None, legend=True):
    probs = np.exp(cdf.iloc[:, 1])
    plt.step(cdf.iloc[:, 0], probs, "k-", label="", clip_on=False, zorder=2)
    ylim = plt.ylim()
    plot_estimate(ylim, expected_value, credible_interval, legend=legend)
    if credible_interval is not None or expected_value is not None:
        _mask_estimate(ylim, cdf.iloc[:, 0], probs, step="pre")
    plt.ylabel("CDF")


def plot_pmf(cdf, expected_value=None, credible_interval=None, legend=True):
    probs = np.exp(cdf.iloc[:, 1])
    probs[1:] = probs[1:] - probs[:-1]
    plt.plot(cdf.iloc[:, 0], probs, "ko", label="", ms=4, clip_on=False, zorder=2)
    ylim = plt.ylim()
    plot_estimate(ylim, expected_value, credible_interval, legend=legend)
    if credible_interval is not None or expected_value is not None:
        _mask_estimate(ylim, cdf.iloc[:, 0], probs)
    plt.ylabel("PMF")


def _mask_estimate(ylim, x, y, step=None):
    plt.fill_between(x, y, ylim[1],  zorder=1, facecolor="white", edgecolor="white", step=step, clip_on=True)
    plt.ylim(ylim)
