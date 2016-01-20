import pandas as pd

from merfishtools.utils import estimate


def read_pmf(path):
    return pd.read_table(path, index_col=[0, 1])


def normalize(*pmfs):
    """Perform a scaling normalization on the given PMFs."""
    if len(pmfs) == 1:
        return pmfs
    expected_values = [estimate(pmf)["ev"] for pmf in pmfs]
    global_median = pd.concat(expected_values).median()
    medians = np.array([ev.median() for ev in expected_values])
    scales = medians / global_median
    for scale, pmf in zip(scales, pmfs):
        pmf["value"] *= scale
