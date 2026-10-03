"""Prespecified modelling engine shared by notebooks 02 and 03.

Everything here was fixed before any model performance was examined:
  - outcome modelled on the logit scale, metrics reported on the fraction scale;
  - elastic net (tuned by inner 5-fold CV) and random forest (fixed settings);
  - 5-fold cross-validation repeated 20 times with seed 2026, identical for every arm;
  - nested linear recalibration and the Nadeau & Bengio (2003) corrected t-test.

The random forest uses RF_TREES = 500 trees (Section 3.5 of the report).
"""
import numpy as np
from scipy.stats import spearmanr, pearsonr, t as tdist
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV, LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RepeatedKFold, KFold, cross_val_predict
from common import logit, expit

SEED = 2026
N_SPLITS = 5
N_REPEATS = 20
RF_TREES = 500


def make_model(kind, seed=0):
    """Return an unfitted model of the given family ('EN' or 'RF')."""
    if kind == "EN":
        # alphas=50 means "a 50-value regularisation path" (scikit-learn >= 1.7 syntax)
        return make_pipeline(
            StandardScaler(),
            ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9, 1.0], alphas=50,
                         cv=KFold(5, shuffle=True, random_state=seed), max_iter=20000))
    if kind == "RF":
        # Fixed, not tuned: tuning more hyperparameters on ~69 samples adds variance.
        # n_jobs=-1 uses all CPU cores; results do not depend on the number of cores.
        return RandomForestRegressor(n_estimators=RF_TREES, min_samples_leaf=3,
                                     max_features=1 / 3, random_state=seed, n_jobs=-1)
    raise ValueError(kind)


def splits(n, n_repeats=None, seed=SEED):
    """The shared train/test splits. Every arm uses the same list, so comparisons are paired."""
    n_repeats = N_REPEATS if n_repeats is None else n_repeats
    return list(RepeatedKFold(n_splits=N_SPLITS, n_repeats=n_repeats,
                              random_state=seed).split(np.zeros(n)))


def fit_predict(kind, Xtr, ytr, Xte_list, recal=False, seed=0):
    """Fit on the logit scale; return back-transformed predictions for each test matrix.

    If recal=True, a linear recalibration map is learned from INNER cross-validated
    predictions on the training data only, and applied to the test predictions. The
    recalibration therefore never sees the samples it corrects.
    """
    mdl = make_model(kind, seed).fit(Xtr, logit(ytr))
    raw = [expit(mdl.predict(X)) for X in Xte_list]
    if not recal:
        return raw, None
    inner = expit(cross_val_predict(make_model(kind, seed), Xtr, logit(ytr),
                                    cv=KFold(5, shuffle=True, random_state=seed + 1)))
    lr = LinearRegression().fit(inner.reshape(-1, 1), ytr)
    rec = [np.clip(lr.predict(p.reshape(-1, 1)), 0.001, 0.999) for p in raw]
    return raw, rec


def pooled_metrics(y, p):
    """Accuracy, agreement and calibration of predictions p against observed y."""
    res = y - p
    slope, intercept = np.polyfit(p, y, 1)      # calibration: regress observed on predicted
    return dict(rmse=float(np.sqrt(np.mean(res ** 2))),
                mae=float(np.mean(np.abs(res))),
                r2=float(1 - np.sum(res ** 2) / np.sum((y - y.mean()) ** 2)),
                spearman=float(spearmanr(y, p)[0]),
                pearson=float(pearsonr(y, p)[0]),
                cal_slope=float(slope), cal_intercept=float(intercept),
                citl=float(np.mean(res)))        # calibration-in-the-large; >0 = underestimation


def nb_corrected_t(diffs, n_train, n_test):
    """Nadeau & Bengio (2003) corrected resampled t-test on per-split differences.

    Overlapping training sets make per-split differences correlated, so the naive
    variance is inflated by (1/J + n_test/n_train) before computing t, CI and p.
    """
    d = np.asarray(diffs)
    J = d.size
    var = np.var(d, ddof=1) * (1 / J + n_test / n_train)
    tstat = d.mean() / np.sqrt(var)
    p = 2 * tdist.sf(abs(tstat), J - 1)
    half = tdist.ppf(0.975, J - 1) * np.sqrt(var)
    return dict(mean=float(d.mean()), lo=float(d.mean() - half), hi=float(d.mean() + half),
                t=float(tstat), p=float(p))
