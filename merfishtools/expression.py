import pandas as pd

from merfishtools.utils import estimate


def read_pmf(path):
    return pd.read_table(path, index_col=[0, 1])


def expected_values(pmf):
    """Calculate matrix of expected values."""
    return (pmf["expr"] * pmf["prob"].exp()).groupby(level=["cell", "feat"]).sum().unstack(0).fillna(0)


def normalize(*pmfs):
    """Perform an upper quartile scaling normalization on the given PMFs."""
    if len(pmfs) == 1:
        return
    # calculate upper quartiles
    quartiles = np.array([expected_values(pmf).unstack().quantile(0.75) for pmf in pmfs])
    mean_quartile = quartiles.mean()
    # calculate scale factors for upper quartile normalization
    scales = [mean_quartile / q for q in quartiles]
    for pmf, scale in zip(pmfs, scales):
        pmf["expr"] *= scale
