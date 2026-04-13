"""
mTBI DTA Meta-analysis — Dataset v15 (all studies verified from full text)
Bivariate pooling approximation via univariate DerSimonian-Laird on logit scale
Date: 2026-04-11
"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# DEFINITIVE DATASET v15 — ALL VERIFIED FROM FULL TEXT
# ============================================================

# --- POOL GFAP+UCH-L1 (7 studies) ---
gfap_studies = [
    # (name, TP, FP, FN, TN, cutoff, platform)
    ("Chayoua 2024",      57, 157,  2,  37, "30/360 i-STAT", "i-STAT POC"),
    ("Bazarian 2021",    115,1061,  5, 720, "30/360 i-STAT", "i-STAT"),
    ("Papa 2022",         23, 245,  0,  81, "67/189 ELISA",  "ELISA Banyan"),
    ("Lagares 2024",     176, 945,  3, 314, "22/327 VIDAS",  "VIDAS bioMérieux"),
    ("Puravet 2026",      50, 700,  0, 260, "35/400 Alinity", "Alinity"),
    ("Milevoj 2025",     107, 529,  5, 181, "35/400 Alinity", "Alinity"),
    ("Legramante 2024",    7,  89,  0,  34, "35/400 Alinity", "Alinity"),
]

# Subgroup: Abbott only (i-STAT + Alinity)
gfap_abbott = [s for s in gfap_studies if "i-STAT" in s[5] or "Alinity" in s[5]]

# --- POOL S100B (4 studies) ---
s100b_studies = [
    # (name, TP, FP, FN, TN, cutoff, window)
    ("Puravet 2026",      48, 720,  2, 240, "0.10", "≤3h"),
    ("Seidenfaden 2021",  32, 366,  0, 168, "0.10", "≤6h"),
    ("Rogan 2023",        15,  81,  1,  36, "0.10", "≤6h"),
    ("Hopman 2023",       69, 355,  5,  66, "0.105","n.s."),
]

# ============================================================
# META-ANALYSIS FUNCTIONS
# ============================================================
def continuity_correction(tp, fp, fn, tn, cc=0.5):
    """Apply continuity correction if any cell is zero"""
    if tp == 0 or fp == 0 or fn == 0 or tn == 0:
        return tp + cc, fp + cc, fn + cc, tn + cc
    return tp, fp, fn, tn

def logit(p):
    return np.log(p / (1 - p))

def inv_logit(x):
    return 1 / (1 + np.exp(-x))

def dl_random_effects(estimates, variances):
    """DerSimonian-Laird random effects on logit scale"""
    w = 1.0 / variances
    w_sum = np.sum(w)
    theta_fe = np.sum(w * estimates) / w_sum
    Q = np.sum(w * (estimates - theta_fe)**2)
    k = len(estimates)
    C = w_sum - np.sum(w**2) / w_sum
    tau2 = max(0, (Q - (k - 1)) / C)
    w_re = 1.0 / (variances + tau2)
    w_re_sum = np.sum(w_re)
    theta_re = np.sum(w_re * estimates) / w_re_sum
    se_re = np.sqrt(1.0 / w_re_sum)
    ci_lo = theta_re - 1.96 * se_re
    ci_hi = theta_re + 1.96 * se_re
    # I-squared
    if Q > 0 and k > 1:
        I2 = max(0, (Q - (k-1)) / Q * 100)
    else:
        I2 = 0
    return {
        'estimate': inv_logit(theta_re),
        'ci_lo': inv_logit(ci_lo),
        'ci_hi': inv_logit(ci_hi),
        'tau2': tau2,
        'I2': I2,
        'Q': Q,
        'k': k
    }

def pool_dta(studies, label=""):
    """Pool sensitivity and specificity separately using DL random effects"""
    cc = 0.5
    se_logits, se_vars = [], []
    sp_logits, sp_vars = [], []
    total_n, total_tp, total_fp, total_fn, total_tn = 0, 0, 0, 0, 0

    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"{'='*70}")
    print(f"{'Study':<22} {'N':>5} {'TP':>5} {'FP':>5} {'FN':>4} {'TN':>5} {'Se%':>6} {'Sp%':>6}")
    print("-"*60)

    for s in studies:
        name, tp, fp, fn, tn = s[0], s[1], s[2], s[3], s[4]
        n = tp + fp + fn + tn
        total_n += n; total_tp += tp; total_fp += fp; total_fn += fn; total_tn += tn
        se_raw = tp / (tp + fn) if (tp+fn) > 0 else 0
        sp_raw = tn / (tn + fp) if (tn+fp) > 0 else 0
        print(f"{name:<22} {n:>5} {tp:>5} {fp:>5} {fn:>4} {tn:>5} {se_raw*100:>5.1f} {sp_raw*100:>5.1f}")

        tpc, fpc, fnc, tnc = continuity_correction(tp, fp, fn, tn, cc)
        se = tpc / (tpc + fnc)
        sp = tnc / (tnc + fpc)
        se_logits.append(logit(se))
        sp_logits.append(logit(sp))
        se_vars.append(1/tpc + 1/fnc)
        sp_vars.append(1/tnc + 1/fpc)

    print("-"*60)
    print(f"{'TOTALE':<22} {total_n:>5} {total_tp:>5} {total_fp:>5} {total_fn:>4} {total_tn:>5}")

    se_logits = np.array(se_logits)
    se_vars = np.array(se_vars)
    sp_logits = np.array(sp_logits)
    sp_vars = np.array(sp_vars)

    se_result = dl_random_effects(se_logits, se_vars)
    sp_result = dl_random_effects(sp_logits, sp_vars)

    # Clinical metrics
    prevalence = (total_tp + total_fn) / total_n
    se_pooled = se_result['estimate']
    sp_pooled = sp_result['estimate']
    npv = (1 - prevalence) * sp_pooled / ((1 - prevalence) * sp_pooled + prevalence * (1 - se_pooled))
    ct_avoided = sp_pooled * (1 - prevalence)
    missed_per_1000 = (1 - se_pooled) * prevalence * 1000

    print(f"\n--- Stime Pooled (DerSimonian-Laird, logit scale) ---")
    print(f"Sensibilità:  {se_result['estimate']*100:.1f}% [{se_result['ci_lo']*100:.1f}–{se_result['ci_hi']*100:.1f}%]")
    print(f"  τ²={se_result['tau2']:.4f}  I²={se_result['I2']:.1f}%  Q={se_result['Q']:.2f} (k={se_result['k']})")
    print(f"Specificità:  {sp_result['estimate']*100:.1f}% [{sp_result['ci_lo']*100:.1f}–{sp_result['ci_hi']*100:.1f}%]")
    print(f"  τ²={sp_result['tau2']:.4f}  I²={sp_result['I2']:.1f}%  Q={sp_result['Q']:.2f} (k={sp_result['k']})")
    print(f"\nPrevalenza CT+: {prevalence*100:.1f}%")
    print(f"NPV: {npv*100:.2f}%")
    print(f"CT evitate: {ct_avoided*100:.1f}% ({ct_avoided*1000:.0f}/1000)")
    print(f"Lesioni mancate: {missed_per_1000:.1f}/1000")

    return {
        'se': se_result, 'sp': sp_result,
        'npv': npv, 'ct_avoided': ct_avoided,
        'missed': missed_per_1000, 'prevalence': prevalence,
        'N': total_n, 'k': len(studies)
    }

# ============================================================
# RUN ANALYSES
# ============================================================
print("=" * 70)
print("  mTBI DTA META-ANALYSIS — DATASET v15 DEFINITIVO")
print("  Tutti i dati 2×2 verificati dai full text PDF")
print("  Data: 11 aprile 2026")
print("=" * 70)

# 1. GFAP+UCH-L1 — Pool primario (7 studi)
gfap_all = pool_dta(gfap_studies, "GFAP+UCH-L1 — Pool primario (7 studi)")

# 2. GFAP+UCH-L1 — Subgroup Abbott only (5 studi)
gfap_abb = pool_dta(gfap_abbott, "GFAP+UCH-L1 — Subgroup Abbott only (5 studi)")

# 3. S100B — Pool primario (4 studi)
s100b_all = pool_dta(s100b_studies, "S100B — Pool primario (4 studi)")

# 4. Head-to-head comparison
print(f"\n{'='*70}")
print(f"  CONFRONTO HEAD-TO-HEAD")
print(f"{'='*70}")
print(f"{'Metrica':<25} {'GFAP+UCH-L1':>15} {'S100B':>15} {'Δ':>10}")
print("-"*65)
for metric, g, s, fmt in [
    ("Sensibilità pooled", gfap_all['se']['estimate']*100, s100b_all['se']['estimate']*100, ".1f"),
    ("Specificità pooled", gfap_all['sp']['estimate']*100, s100b_all['sp']['estimate']*100, ".1f"),
    ("NPV", gfap_all['npv']*100, s100b_all['npv']*100, ".2f"),
    ("CT evitate/1000", gfap_all['ct_avoided']*1000, s100b_all['ct_avoided']*1000, ".0f"),
    ("Missed/1000", gfap_all['missed'], s100b_all['missed'], ".1f"),
]:
    delta = g - s
    print(f"{metric:<25} {g:>15{fmt}} {s:>15{fmt}} {delta:>+10{fmt}}")
print(f"{'N totale':<25} {gfap_all['N']:>15} {s100b_all['N']:>15}")
print(f"{'k studi':<25} {gfap_all['k']:>15} {s100b_all['k']:>15}")
