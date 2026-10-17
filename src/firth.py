"""Firth-penalised logistic regression for aggregated binomial data, and tests for a trend in log RORs.

The modified score of Firth (1993) is U* = X'(y - n*p + h*(1/2 - p)), where h is the diagonal of the hat
matrix under weights w = n*p*(1-p). Wald 95% CIs are computed from (X'WX)^-1. The penalty keeps the
estimates finite when a cell is empty.
"""
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import chi2, norm


def firth_binom(X, y, n, maxiter=100, tol=1e-8):
    """X: design matrix (with intercept), y: events per stratum, n: reports per stratum.
    Returns (coefficients, standard errors)."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n = np.asarray(n, float)
    beta = np.zeros(X.shape[1])
    beta[0] = np.log((y.sum() + 0.5) / (n.sum() - y.sum() + 0.5))
    for _ in range(maxiter):
        p = expit(X @ beta)
        w = n * p * (1 - p)
        inv = np.linalg.pinv(X.T @ (w[:, None] * X))
        h = w * np.einsum("ij,jk,ik->i", X, inv, X)
        step = inv @ (X.T @ (y - n * p + h * (0.5 - p)))
        mx = np.abs(step).max()
        if mx > 5:  # damp very large steps
            step *= 5 / mx
        beta += step
        if np.abs(step).max() < tol:
            break
    p = expit(X @ beta)
    inv = np.linalg.pinv(X.T @ ((n * p * (1 - p))[:, None] * X))
    return beta, np.sqrt(np.diag(inv))


def fit_ror(data, y, expo, covars, ref):
    """data: one row per report. Adjusted Firth ROR for each exposure column (0/1).
    covars: categorical columns (dummy-coded, reference level from `ref`) or 0/1 columns.
    Returns {exposure: (ROR, lower, upper, logROR, SE)}."""
    cols = expo + covars
    g = data.groupby(cols, observed=True)[y].agg(["sum", "count"]).reset_index()
    cat_cols = [c for c in covars if c in ref]
    X = pd.get_dummies(g[cols], columns=cat_cols)
    X = X.drop(columns=[f"{c}_{v}" for c, v in ref.items() if f"{c}_{v}" in X.columns]).astype(float)
    X = X.loc[:, X.std() > 0]  # drop columns that are constant in this subset
    X.insert(0, "const", 1.0)
    beta, se = firth_binom(X.values, g["sum"].values, g["count"].values)
    idx = {c: i for i, c in enumerate(X.columns)}
    out = {}
    for e in expo:
        b, s = beta[idx[e]], se[idx[e]]
        out[e] = (np.exp(b), np.exp(b - 1.96 * s), np.exp(b + 1.96 * s), b, s)
    return out


def gradient_test(logrors, ses):
    """Inverse-variance weighted regression of log RORs on the tier rank (0, 1, 2, ...), plus Cochran's Q.
    Tiers are treated as independent (reports with events in more than one tier are rare)."""
    b = np.asarray(logrors, float)
    s = np.asarray(ses, float)
    w = 1 / s**2
    k = np.arange(len(b), dtype=float)
    kbar = (w * k).sum() / w.sum()
    slope = (w * (k - kbar) * b).sum() / (w * (k - kbar) ** 2).sum()
    se_slope = np.sqrt(1 / (w * (k - kbar) ** 2).sum())
    z = slope / se_slope
    pooled = (w * b).sum() / w.sum()
    q = (w * (b - pooled) ** 2).sum()
    return dict(slope=slope, se=se_slope, ratio_per_step=np.exp(slope), p_trend=2 * norm.sf(abs(z)),
                Q=q, p_het=chi2.sf(q, len(b) - 1))


def ratio_test(b1, s1, b2, s2):
    """Ratio of two RORs (ROR1/ROR2) with 95% CI and p value (z test, independence assumed)."""
    d = b1 - b2
    se = np.sqrt(s1**2 + s2**2)
    return np.exp(d), np.exp(d - 1.96 * se), np.exp(d + 1.96 * se), 2 * norm.sf(abs(d / se))
