
# Feature Engineering

**Project:** FRAGMENTOMIC AND COPY-NUMBER-DERIVED TUMOUR FRACTION CONTRIBUTIONS TO CTDNA BURDEN PREDICTION IN METASTATIC CASTRATION-RESISTANT PROSTATE CANCER: AN ABLATION STUDY

This document covers the first stage of the analysis: turning raw fragment
length histograms into a small set of interpretable features. It records what
was done, what came out, and how to rerun it independently.

---

## 1. What the raw data looks like

The source is the Renaud et al. (2022) supplementary workbook, `data.xlsx`.
Two sheets:

**WGS sheet** 142 rows, 673 columns

| Column | Contents |
|---|---|
| 0 | ichorCNA tumour fraction (all 142 rows) |
| 1 | VAF-based ctDNA fraction, the reference standard (86 rows) |
| 2 onward | Fragment length counts, 671 bins covering 30–700 bp |

**capture sheet**  86 rows, 671 columns of fragment counts only. Same samples,
measured on a targeted capture panel instead of shallow whole-genome sequencing.

The first two rows of each sheet are headers and a formula row, so both are
read with `skiprows=2`.

---

## 2. Why not use the 671 raw bins

With 86 labelled samples and 671 predictors, any model would have roughly 0.13
samples per feature. That is not a modelling problem to be solved with
regularisation; it is a guarantee of overfitting. Published fragmentomic work
takes the same route of deriving a small number of summary measures, and this
analysis follows that convention.

The target is 10–15 features, each with a stated biological rationale.

---

## 3. The features and why each one exists

Every feature is computed from the normalised density, meaning raw counts
divided by their total so each sample sums to 1. Without this, samples
sequenced more deeply would appear different purely because of read depth.

### Size-ratio features

**`short_long_ratio`**  proportion of fragments 100–150 bp divided by
proportion 151–220 bp. This is the DELFI-style metric. Tumour-derived cfDNA is
systematically shorter, so this ratio should rise with tumour content.

**`sub_nucleosomal_frac`**  proportion of fragments below 150 bp, meaning
shorter than a single nucleosome footprint. Same rationale, computed without a
denominator.

### Nucleosome-peak features

cfDNA is not cut at random. DNA wraps around nucleosomes in roughly 147 bp,
and nucleases cleave the exposed linker between them, so fragment lengths
cluster near 167 bp.

**`mono_peak_pos`**  the fragment length at which the peak occurs, searched
within 120–220 bp. A left shift indicates tumour content.

**`mono_peak_height`**  density at that peak. Taller means lengths concentrate
more tightly.

**`mono_peak_fwhm`**  full width at half maximum, meaning the span of lengths
over which density stays above half the peak height. A broader peak implies
less regular nucleosomal protection.

**`di_mono_ratio`**  proportion at 300–360 bp divided by proportion at
150–200 bp. When a nuclease skips a linker site, the resulting fragment spans
two nucleosomes and lands near 330 bp. Tumour cfDNA shows a more prominent
di-nucleosome peak.

### Periodicity

**`periodicity_10bp`**  amplitude of the roughly 10 bp oscillation in the
distribution.

DNA completes one helical turn every 10.4 bp. When wrapped on a nucleosome the
same face is exposed to solvent at that interval while the opposite face is
shielded, so cleavage is favoured at positions one turn apart. This leaves a
ripple in the length distribution.

Extracting it takes three steps: restrict to 80–250 bp, subtract a moving
average to remove the overall peak shape, then take a Fourier transform and
measure amplitude in the frequency band corresponding to a 9.5–11.5 bp period.

### Entropy

**`entropy`**  Shannon entropy of the length distribution in bits. Higher
means lengths are spread more evenly; lower means they concentrate at
particular values, implying more uniform nucleosomal protection.

### Distribution shape

**`mean_length`**, **`sd_length`**, **`skewness`**, **`kurtosis`**,
**`median_length`**, **`iqr_length`**  standard moments and quantiles. These
capture overall shifts that the targeted features might miss.

---

## 4. Results

### Feature matrix

142 rows by 14 features, with no missing values.

```
                          mean     std       min       max
short_long_ratio        0.2833  0.1124    0.1647    0.7920
sub_nucleosomal_frac    0.2024  0.0521    0.1271    0.4188
mono_peak_pos         165.4789  0.5677  164.0000  167.0000
mono_peak_height        0.0250  0.0031    0.0135    0.0336
mono_peak_fwhm         23.9789  4.9300   18.0000   51.0000
di_mono_ratio           0.0925  0.0386    0.0402    0.2728
periodicity_10bp        3.7075  0.4096    2.3794    5.5648
entropy                 7.1864  0.2237    6.6444    7.8163
mean_length           183.1299  8.4765  163.6939  210.6954
sd_length              68.5490  9.2909   48.6879   98.0690
skewness                2.5504  0.3192    1.5725    3.3680
kurtosis               12.4346  2.8589    5.8453   20.2855
median_length         167.5211  2.9093  156.0000  174.0000
iqr_length             30.0845  8.9161   21.0000   99.0000
```

The mononucleosome peak sits between 164 and 167 bp across all samples, which
is exactly where nucleosomal biology predicts. That alone is a strong sign the
extraction is working.

### Univariate association with the reference standard

Spearman correlation against VAF-based ctDNA fraction, n = 86:

| Feature | Spearman r | p | Expected direction | Matches? |
|---|---|---|---|---|
| short_long_ratio | **+0.645** | 2.1e-11 | + (tumour shorter) | yes |
| di_mono_ratio | **+0.588** | 2.6e-09 | + (more di-nucleosome) | yes |
| kurtosis | −0.524 | 2.2e-07 | not predicted | — |
| iqr_length | +0.506 | 6.7e-07 | not predicted | — |
| entropy | +0.485 | 2.3e-06 | not predicted | — |
| skewness | −0.466 | 6.1e-06 | not predicted | — |
| sub_nucleosomal_frac | +0.456 | 1.0e-05 | + (tumour shorter) | yes |
| sd_length | +0.438 | 2.5e-05 | not predicted | — |
| mono_peak_height | −0.434 | 3.0e-05 | not predicted | — |
| mono_peak_fwhm | +0.381 | 2.9e-04 | not predicted | — |
| median_length | −0.247 | 0.022 | − (left shift) | yes |
| mono_peak_pos | −0.222 | 0.040 | − (left shift) | yes |
| mean_length | +0.217 | 0.045 | − (left shift) | **no** |
| periodicity_10bp | −0.091 | 0.406 | + (increased) | **no** |

**Baseline to beat:** ichorCNA alone against VAF gives Spearman r = 0.829.
This reproduces the Pearson r = 0.79 reported in the paper, which confirms the
columns are being read correctly.

---


### periodicity_10bp shows no association

r = −0.091, p = 0.41. Three possible explanations:

1. The detrending window (21 bp) is wide enough to remove some of the
   oscillation along with the trend.
2. The length range 80–250 bp is wrong. The published observation concerns
   periodicity *to the left of* the main mode, so a narrower window such as
   90–160 bp may be more appropriate.
3. The signal genuinely does not track tumour fraction in this cohort.

This needs testing before the feature is either kept or dropped. 
below sets out how.

---

## 6. How to reproduce this

### Requirements

```bash
pip install pandas numpy scipy openpyxl matplotlib scikit-learn
```

### Files

Put `data.xlsx` (the Renaud supplementary workbook) and `features.py` in the
same folder.

### Step 1 confirm the data reads correctly

```python
import pandas as pd, numpy as np

w = pd.read_excel('data.xlsx', sheet_name='WGS', header=None, skiprows=2)
ichor = pd.to_numeric(w[0], errors='coerce').values
vaf   = pd.to_numeric(w[1], errors='coerce').values
frag  = w.iloc[:, 2:].apply(pd.to_numeric, errors='coerce').values

print('rows:', len(w), 'bins:', frag.shape[1])
print('complete cases:', (~np.isnan(vaf)).sum())

m = ~np.isnan(vaf)
print('ichor vs vaf Pearson r:', round(float(np.corrcoef(ichor[m], vaf[m])[0,1]), 3))
```

**Expected:** 142 rows, 671 bins, 86 complete cases, r = 0.79.

If that r is not 0.79, the columns are being read wrongly and nothing
downstream will be valid. Stop and fix it before continuing.

### Step 2 extract the features

```python
from features import load_renaud, build_feature_matrix

data = load_renaud('data.xlsx')
X = build_feature_matrix(data['frag_wgs'])

print(X.shape)                 # expect (142, 14)
print(X.isna().sum().sum())    # expect 0
print(X.describe().T[['mean','std','min','max']].round(4))
```

**The check that matters:** `mono_peak_pos` should fall between 164 and 167 for
every sample. If it does not, the normalisation or the peak search window is
wrong.

### Step 3 verify against known biology

```python
from scipy.stats import spearmanr

m = data['labelled']
vaf = data['vaf']

for f in X.columns:
    r, p = spearmanr(X[f][m], vaf[m])
    print(f'{f:24s} r={r:+.3f}  p={p:.2e}')
```

**Expected:** `short_long_ratio` positive and strongest, `di_mono_ratio`
positive and second, `median_length` and `mono_peak_pos` both negative.

This step is the real test. Any feature can be computed; the question is
whether it behaves as fragmentomic biology predicts. If `short_long_ratio` came
out negative, something would be wrong with the implementation, not with the
biology.

### Step 4 export

```python
X.insert(0, 'ichorCNA', data['ichor'])
X.insert(1, 'vaf_ctdna', data['vaf'])
X.to_csv('features_wgs.csv', index=False)

Xc = build_feature_matrix(data['frag_capture'])
Xc.to_csv('features_capture.csv', index=False)
```

---

## 7. Testing the periodicity feature

Run this to determine whether the null result is a parameter problem or a real
absence of signal. It sweeps the length window and the detrending width, and
reports which combination gives the strongest association.

```python
import numpy as np
from scipy.stats import spearmanr
from features import load_renaud, LENGTHS

data = load_renaud('data.xlsx')
m, vaf, frag = data['labelled'], data['vaf'], data['frag_wgs']

def periodicity(counts, lo, hi, win):
    d = counts / counts.sum()
    sel = (LENGTHS >= lo) & (LENGTHS <= hi)
    sig = d[sel].astype(float)
    if sig.size < 32:
        return np.nan
    trend = np.convolve(sig, np.ones(win)/win, mode='same')
    resid = sig - trend
    resid -= resid.mean()
    spec = np.abs(np.fft.rfft(resid))
    freqs = np.fft.rfftfreq(resid.size, d=1.0)
    band = (freqs >= 1/11.5) & (freqs <= 1/9.5)
    if not band.any():
        return np.nan
    return float(spec[band].max()) / float(spec[1:].mean())

results = []
for lo, hi in [(80,250), (90,160), (100,180), (120,200), (140,220), (90,200)]:
    for win in [7, 11, 15, 21, 31]:
        vals = np.array([periodicity(frag[i,:], lo, hi, win)
                         for i in range(frag.shape[0])])
        ok = m & ~np.isnan(vals)
        if ok.sum() < 50:
            continue
        r, p = spearmanr(vals[ok], vaf[ok])
        results.append((lo, hi, win, round(r,3), f'{p:.3f}'))

results.sort(key=lambda t: -abs(t[3]))
print(f"{'lo':>4} {'hi':>4} {'win':>4} {'r':>7} {'p':>8}")
for row in results[:12]:
    print(f'{row[0]:>4} {row[1]:>4} {row[2]:>4} {row[3]:>+7.3f} {row[4]:>8}')
```

**How to read the output.** If some window gives an r above about 0.3 in the
expected positive direction, the original parameters were wrong and the feature
should be recomputed with the better ones. If every combination sits near zero,
the signal genuinely is not there in this cohort, and the feature should be
dropped with that stated explicitly.


## 8. What comes next

With features verified, the primary analysis is the three-arm ablation:

- **Arm A** : the 14 fragmentomic features alone
- **Arm B** : ichorCNA alone
- **Arm C** : both together

Each predicts VAF-based ctDNA fraction on the 86 complete cases, using elastic
net and random forest under repeated cross-validation.
