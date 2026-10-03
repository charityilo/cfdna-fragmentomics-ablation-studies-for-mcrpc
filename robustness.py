"""Robustness and explainability analyses.

1. Robustness to the choice of learning algorithm and to hyperparameter tuning: XGBoost and
   support vector regression (SVR), each fitted both with tuning by nested cross-validation
   inside every outer training fold and with its library's default settings.
2. Grouped permutation importance on held-out test folds, with features grouped by the four
   aspects of the length distribution defined in the report (Section 3.5), plus ichorCNA.
Everything uses the same 100 outer splits, outcome transformation and metrics as the ablation
in notebook 02, so the results are directly comparable.
"""
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV, KFold
from xgboost import XGBRegressor
import engine
from common import logit, expit

# Tuning grids for XGBoost and SVR
XGB_GRID = {"max_depth": [1, 2, 3], "learning_rate": [0.03, 0.1], "n_estimators": [100, 300]}
SVR_GRID = {"svr__C": [0.1, 1, 10, 100], "svr__gamma": ["scale", 0.01, 0.1], "svr__epsilon": [0.05, 0.2]}

# Feature groups: the four aspects of the fragment-length distribution described in the report
GROUPS = {
    "Short-fragment enrichment": ["short_long_ratio", "sub_nucleosomal_frac"],
    "Mononucleosomal peak shape": ["mono_peak_pos", "mono_peak_height", "mono_peak_fwhm"],
    "Dinucleosomal content": ["di_mono_ratio"],
    "Overall distribution shape": ["periodicity_10bp", "entropy", "mean_length", "sd_length",
                                   "skewness", "kurtosis", "median_length", "iqr_length"],
}


def make_extra_model(kind, seed=0):
    """A tuned XGBoost or SVR: GridSearchCV with an inner five-fold split of the training fold."""
    inner = KFold(5, shuffle=True, random_state=seed)
    if kind == "XGB":
        base = XGBRegressor(subsample=0.8, colsample_bytree=0.8, objective="reg:squarederror",
                            random_state=seed, n_jobs=1, verbosity=0)
        return GridSearchCV(base, XGB_GRID, cv=inner, scoring="neg_mean_squared_error", n_jobs=1)
    if kind == "SVR":
        base = make_pipeline(StandardScaler(), SVR(kernel="rbf"))
        return GridSearchCV(base, SVR_GRID, cv=inner, scoring="neg_mean_squared_error", n_jobs=1)
    # Untuned versions: the libraries' default settings
    if kind == "XGB_default":
        return XGBRegressor(objective="reg:squarederror", random_state=seed, n_jobs=1, verbosity=0)
    if kind == "SVR_default":
        return make_pipeline(StandardScaler(), SVR(kernel="rbf"))
    raise ValueError(kind)


def fit(kind, Xtr, ytr, seed):
    """Fit on the logit scale. EN and RF are the primary models from engine.py."""
    mdl = make_extra_model(kind, seed) if kind in ("XGB", "SVR", "XGB_default", "SVR_default") else engine.make_model(kind, seed)
    return mdl.fit(Xtr, logit(ytr))


def rmse(y, p):
    return float(np.sqrt(np.mean((y - p) ** 2)))


def grouped_permutation_importance(mdl, Xte, yte, groups_idx, n_perm=5, seed=0):
    """Increase in test-fold RMSE when all columns of a group are shuffled together.

    Shuffling a whole group at once avoids splitting credit between near-duplicate features.
    """
    rng = np.random.default_rng(seed)
    base = rmse(yte, expit(mdl.predict(Xte)))
    out = {}
    for name, cols in groups_idx.items():
        inc = []
        for _ in range(n_perm):
            Xp = Xte.copy()
            Xp[:, cols] = Xte[rng.permutation(len(Xte))][:, cols]
            inc.append(rmse(yte, expit(mdl.predict(Xp))) - base)
        out[name] = float(np.mean(inc))
    return out


def run(kind, arm, X, y, S, groups_idx=None):
    """Fit one model family to one arm over all splits; optionally compute importance."""
    raw = np.full((engine.N_REPEATS, len(y)), np.nan)
    fold_rmse, params, importance, coefs = [], [], [], []
    for s, (tr, te) in enumerate(S):
        mdl = fit(kind, X[tr], y[tr], seed=s)
        p = expit(mdl.predict(X[te]))
        raw[s // engine.N_SPLITS, te] = p
        fold_rmse.append(rmse(y[te], p))
        if hasattr(mdl, "best_params_"):
            params.append(mdl.best_params_)
        if kind == "EN":
            coefs.append(mdl[-1].coef_.copy())
        if groups_idx is not None:
            importance.append(grouped_permutation_importance(mdl, X[te], y[te], groups_idx, seed=s))
    return dict(kind=kind, arm=arm, raw=raw, fold_rmse=np.array(fold_rmse), params=params,
                importance=importance, coefs=coefs)


def holm(pvals):
    """Holm step-down adjustment for a small family of p-values."""
    order = np.argsort(pvals); m = len(pvals); adj = np.empty(m); running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * pvals[i])); adj[i] = running
    return adj
