import pandas as pd


def estimate(pmf):
    values = lambda pmf: pmf.iloc[:, 0]
    probs = lambda pmf: pmf.iloc[:, 1]

    def agg(pmf):
        ev = values(pmf) * probs(pmf).exp().sum()
        var = ((values(feat_pmf) - ev) ** 2 * probs(feat_pmf)).sum()
        std = var.sqrt()
        return pd.DataFrame({"ev": ev, "var": var, "std": std})

    return pmf.group_by(level=["cell", "feat"]).agg(agg)
