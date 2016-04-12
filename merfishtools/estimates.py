import pandas as pd
import matplotlib.pyplot as plt


def read_exp_estimates(path_or_buffer):
    return pd.read_table(path_or_buffer, index_col=[0, 1])


def read_diffexp_estimates(path_or_buffer):
    return pd.read_table(path_or_buffer, index_col=0)


def plot_estimate(ylim, expected_value=None, credible_interval=None, legend=True):
    if credible_interval is not None:
        ci_lower, ci_upper = credible_interval
        plt.fill([ci_lower, ci_upper, ci_upper, ci_lower], [0, 0, ylim[1], ylim[1]], "red", lw=0, label="95% credible interval", alpha=0.5)
    if expected_value is not None:
        plt.vlines([expected_value], *ylim, colors="red", linestyles="-", label="expected value")
    if legend:
        plt.legend(loc="best")
    plt.ylim(ylim)
