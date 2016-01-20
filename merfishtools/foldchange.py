import pandas as pd


def read_pmf(path):
    return pd.read_table(path, index_col=0)
