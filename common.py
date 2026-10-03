"""Data loading and feature extraction shared by all three notebooks.

Keeping these functions in one module (rather than copying them into each notebook)
guarantees that every analysis uses exactly the same definitions.
"""
import hashlib
import numpy as np
import pandas as pd

# Fragment lengths 30-700 bp, one column per base pair (671 bins). Row 2 of each sheet
# holds these lengths, which is why both sheets are read with skiprows=2.
LENGTHS = np.arange(30, 701)

# Checksum of Supplementary File 1, version 2 (Renaud et al., 2022). If the file on disk
# does not match, it is not the file analysed in the report.
EXPECTED_SHA256 = "f0693eee0b6b2c0ddcdcfd609ba9ab087828db73313a135b05afb4859c4698e6"


def sha256(path):
    """Return the SHA-256 checksum of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load(path):
    """Read both sheets of the supplementary workbook.

    Returns a dict with the ichorCNA estimate, the VAF-based ctDNA fraction (NaN where the
    file says "NA"), the raw sWGS and capture fragment counts, and a mask of the samples
    that have a VAF-based reference value.
    """
    w = pd.read_excel(path, sheet_name="WGS", header=None, skiprows=2)
    c = pd.read_excel(path, sheet_name="capture", header=None, skiprows=2)
    ichor = pd.to_numeric(w[0], errors="coerce").to_numpy(float)
    vaf = pd.to_numeric(w[1], errors="coerce").to_numpy(float)
    frag_wgs = w.iloc[:, 2:].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    frag_cap = c.apply(pd.to_numeric, errors="coerce").to_numpy(float)
    return dict(ichor=ichor, vaf=vaf, frag_wgs=frag_wgs, frag_cap=frag_cap,
                labelled=~np.isnan(vaf))


def _frac(d, lo, hi):
    """Proportion of fragments with length in [lo, hi] bp."""
    sel = (LENGTHS >= lo) & (LENGTHS <= hi)
    return d[sel].sum()


def _weighted_quantile(d, q):
    """Fragment length at which the cumulative density first reaches q."""
    return float(LENGTHS[np.searchsorted(np.cumsum(d), q)])


def periodicity(d, lo=80, hi=250, win=21):
    """Strength of the ~10-bp oscillation in the length distribution.

    The density between lo and hi bp is detrended with a moving average, so only the
    fine oscillation remains; the largest Fourier amplitude at periods of 9.5-11.5 bp is
    then divided by the mean amplitude over all non-zero frequencies. The window and
    detrending span follow the feature specification and were NOT tuned on the outcome.
    """
    sel = (LENGTHS >= lo) & (LENGTHS <= hi)
    sig = d[sel].astype(float)
    trend = np.convolve(sig, np.ones(win) / win, mode="same")
    resid = sig - trend
    resid -= resid.mean()
    spec = np.abs(np.fft.rfft(resid))
    freqs = np.fft.rfftfreq(resid.size, d=1.0)
    band = (freqs >= 1 / 11.5) & (freqs <= 1 / 9.5)
    return float(spec[band].max()) / float(spec[1:].mean())


def features_one(counts):
    """The 14 fragmentomic features (Table 3.2) for one sample's fragment counts."""
    d = counts / counts.sum()          # density: removes sequencing depth as a signal
    f = {}
    f["short_long_ratio"] = _frac(d, 100, 150) / _frac(d, 151, 220)   # DELFI-style ratio
    f["sub_nucleosomal_frac"] = _frac(d, 30, 149)
    # Mononucleosomal peak, searched between 120 and 220 bp
    win = (LENGTHS >= 120) & (LENGTHS <= 220)
    dw, lw = d[win], LENGTHS[win]
    ip = int(np.argmax(dw))
    f["mono_peak_pos"] = float(lw[ip])
    f["mono_peak_height"] = float(dw[ip])
    half = dw[ip] / 2               # width: contiguous run of bins at >= half the peak height
    left = ip
    while left - 1 >= 0 and dw[left - 1] >= half:
        left -= 1
    right = ip
    while right + 1 < len(dw) and dw[right + 1] >= half:
        right += 1
    f["mono_peak_fwhm"] = float(lw[right] - lw[left] + 1)
    f["di_mono_ratio"] = _frac(d, 300, 360) / _frac(d, 150, 200)
    f["periodicity_10bp"] = periodicity(d)
    nz = d[d > 0]
    f["entropy"] = float(-(nz * np.log2(nz)).sum())
    # Moments of the length distribution
    mu = float((LENGTHS * d).sum())
    f["mean_length"] = mu
    f["sd_length"] = float(np.sqrt(((LENGTHS - mu) ** 2 * d).sum()))
    sd = f["sd_length"]
    f["skewness"] = float((((LENGTHS - mu) / sd) ** 3 * d).sum())
    f["kurtosis"] = float((((LENGTHS - mu) / sd) ** 4 * d).sum()) - 3.0   # EXCESS kurtosis
    f["median_length"] = _weighted_quantile(d, 0.5)
    f["iqr_length"] = _weighted_quantile(d, 0.75) - _weighted_quantile(d, 0.25)
    return f


def feature_matrix(frag):
    """Apply features_one to every row of a count matrix."""
    return pd.DataFrame([features_one(r) for r in frag])


def logit(p):
    """Logit with clipping, so fractions of exactly 0 or 1 cannot produce infinities."""
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return np.log(p / (1 - p))


def expit(x):
    """Inverse logit: maps model outputs back to the fraction scale."""
    return 1 / (1 + np.exp(-x))
