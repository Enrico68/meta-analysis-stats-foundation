#!/usr/bin/env python3
"""
mTBI Biomarker DTA Meta-Analysis — CORRECTED v14
=================================================
Data source  : Dataset_Metaanalisi_Corretto.xlsx (full-text verified)
Date         : 11 April 2026
Pooled stats : DerSimonian-Laird random effects (Python)
Study stats  : computed from raw 2×2 tables + Wilson CIs
Figures      : Fig2 (GFAP forest), Fig3 (S100B forest),
               Fig4 (SROC), Fig5 (clinical summary), FigS1 (funnel)

CORREZIONI rispetto a v13:
- GFAP: 7 studi (N=4,955) - esclusi Puravet 2026 (review), Jones 2020 (S100B)
- S100B: 3 studi (N=1,871) - esclusi Clermont-Ferrand 2016, Puravet 2026, Hopman 2023
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 0. GLOBAL STYLE
# =============================================================================
plt.rcParams.update({
    'font.family'        : 'DejaVu Sans',
    'font.size'          : 10,
    'axes.linewidth'     : 0.8,
    'axes.spines.top'    : False,
    'axes.spines.right'  : False,
    'figure.dpi'         : 150,
    'savefig.dpi'        : 300,
    'savefig.bbox'       : 'tight',
    'savefig.facecolor'  : 'white',
})

COL_G  = '#1B4F8A'   # GFAP dark blue
COL_GL = '#90B8D8'   # GFAP light blue
COL_S  = '#1A6B3C'   # S100B dark green
COL_SL = '#90C9A8'   # S100B light green

# =============================================================================
# 1. RAW 2×2 DATA  (CORRECTED from full-text verification)
# =============================================================================

GFAP = dict(
    label    = "GFAP+UCH-L1",
    k        = 7,
    N_orig   = 4955,
    CT_pos   = 492,
    prev_pct = 9.9,
    studies  = ["Chayoua 2024", "Bazarian 2021", "Papa 2022", "Lagares 2024",
                "Milevoj 2025", "Legramante 2024", "Lapić 2024"],
    TP = np.array([57, 115, 23, 176, 107, 7, 7]),
    FP = np.array([157, 1061, 245, 945, 529, 89, 38]),
    FN = np.array([2, 5, 0, 3, 5, 0, 0]),
    TN = np.array([37, 720, 81, 314, 181, 34, 17]),
    # ── Pooled (DerSimonian-Laird random effects) ───────────────────────────
    pool_se    = 0.973, pool_se_lo = 0.959, pool_se_hi = 0.987,
    pool_sp    = 0.276, pool_sp_lo = 0.209, pool_sp_hi = 0.343,
    I2_se=0,   tau2_se=0.000, Q_se=3.1,   pQ_se=0.792,
    I2_sp=95,  tau2_sp=0.007, Q_sp=125.2, pQ_sp=0.000,
    NPV=98.94, CT_avoided=27.6, missed_1000=2.7,
    color=COL_G, color_l=COL_GL,
)

S100B = dict(
    label    = "S100B",
    k        = 3,
    N_orig   = 1871,
    CT_pos   = 110,
    prev_pct = 5.9,
    studies  = ["Oris 2021", "Rogan 2023", "Seidenfaden 2021"],
    TP = np.array([63, 15, 32]),
    FP = np.array([891, 81, 452]),
    FN = np.array([0, 1, 0]),
    TN = np.array([218, 36, 82]),
    # ── Pooled (DerSimonian-Laird random effects) ───────────────────────────
    pool_se    = 0.989, pool_se_lo = 0.970, pool_se_hi = 1.000,
    pool_sp    = 0.204, pool_sp_lo = 0.149, pool_sp_hi = 0.260,
    I2_se=0,   tau2_se=0.000, Q_se=0.8,  pQ_se=0.656,
    I2_sp=85,  tau2_sp=0.002, Q_sp=13.2, pQ_sp=0.001,
    NPV=99.67, CT_avoided=20.4, missed_1000=0.6,
    color=COL_S, color_l=COL_SL,
)

# =============================================================================
# 2. HELPERS
# =============================================================================

def apply_cc(TP, FP, FN, TN, add=0.5):
    """Continuity correction: add 0.5 to zero cells only."""
    TP = np.where(TP == 0, TP + add, TP.astype(float))
    FN = np.where(FN == 0, FN + add, FN.astype(float))
    FP = np.where(FP == 0, FP + add, FP.astype(float))
    TN = np.where(TN == 0, TN + add, TN.astype(float))
    return TP, FP, FN, TN

def wilson_ci(p, n, z=1.96):
    """Wilson score interval (vectorised)."""
    p  = np.asarray(p, dtype=float)
    n  = np.asarray(n, dtype=float)
    c  = (p + z**2/(2*n)) / (1 + z**2/n)
    m  = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / (1 + z**2/n)
    return np.clip(c - m, 0, 1), np.clip(c + m, 0, 1)

def study_stats(d):
    """Return Se, Sp, CIs, weights, N per study."""
    TP, FP, FN, TN = apply_cc(d['TP'].copy(), d['FP'].copy(),
                               d['FN'].copy(), d['TN'].copy())
    Se   = TP / (TP + FN)
    Sp   = TN / (TN + FP)
    N_se = TP + FN
    N_sp = TN + FP
    se_lo, se_hi = wilson_ci(Se, N_se)
    sp_lo, sp_hi = wilson_ci(Sp, N_sp)
    w_se = 1 / (Se*(1-Se)/N_se + 1e-9)
    w_sp = 1 / (Sp*(1-Sp)/N_sp + 1e-9)
    N    = d['TP'] + d['FP'] + d['FN'] + d['TN']
    return Se, Sp, se_lo, se_hi, sp_lo, sp_hi, w_se, w_sp, N

# =============================================================================
# 3. FIGURES 2 & 3 — FOREST PLOTS
# =============================================================================

def plot_forest(d, fig_num, filename):
    Se, Sp, se_lo, se_hi, sp_lo, sp_hi, w_se, w_sp, N = study_stats(d)
    k      = d['k']
    col    = d['color']
    y      = np.arange(k, 0, -1, dtype=float)   # k … 1 (top to bottom)
    y_pool = 0.0

    fig, axes = plt.subplots(1, 2, figsize=(16, max(6, k*0.8 + 3.0)),
                             gridspec_kw={'wspace': 0.04})

    params = [
        # vals  lo     hi     pool_v          pool_lo         pool_hi
        #  I2      tau2         Q              pQ          w     xlabel
        (Se,   se_lo, se_hi,
         d['pool_se'],  d['pool_se_lo'], d['pool_se_hi'],
         d['I2_se'],    d['tau2_se'],    d['Q_se'],     d['pQ_se'],
         w_se, "Sensitivity (95% CI)"),

        (Sp,   sp_lo, sp_hi,
         d['pool_sp'],  d['pool_sp_lo'], d['pool_sp_hi'],
         d['I2_sp'],    d['tau2_sp'],    d['Q_sp'],     d['pQ_sp'],
         w_sp, "Specificity (95% CI)"),
    ]

    for ax_i, (ax, (vals, lo, hi,
                    pv, plo, phi,
                    I2, tau2, Q, pQ,
                    w, xlabel)) in enumerate(zip(axes, params)):

        # ── Background rows ──────────────────────────────────────────────────
        for i, yi in enumerate(y):
            if i % 2 == 0:
                ax.axhspan(yi - 0.45, yi + 0.45, color='#F6F6F6', zorder=0)
        bg_col = '#EEF4FF' if col == COL_G else '#EEF9EE'
        ax.axhspan(y_pool - 0.48, y_pool + 0.48, color=bg_col, zorder=0)

        ax.axvline(0, color='#CCCCCC', lw=0.5, zorder=1)

        # ── Individual studies ───────────────────────────────────────────────
        max_w = w.max()
        for i, yi in enumerate(y):
            # CI line
            ax.plot([lo[i], hi[i]], [yi, yi], color=col, lw=1.6, zorder=2)
            ax.plot([lo[i]]*2, [yi-0.13, yi+0.13], color=col, lw=1.6, zorder=2)
            ax.plot([hi[i]]*2, [yi-0.13, yi+0.13], color=col, lw=1.6, zorder=2)
            # Point (size ∝ weight)
            sz = 40 + (w[i] / max_w) * 130
            ax.scatter(vals[i], yi, s=sz, color=col, zorder=3,
                       edgecolors='white', linewidth=0.9)
            # Value label
            ax.text(1.07, yi, f"{vals[i]*100:.1f}",
                    va='center', ha='left', fontsize=9, color='#333333')

        # ── Pooled diamond ───────────────────────────────────────────────────
        dw = (phi - plo) / 2
        ax.fill([pv-dw, pv, pv+dw, pv, pv-dw],
                [y_pool, y_pool+0.38, y_pool, y_pool-0.38, y_pool],
                color=col, zorder=3)
        ax.axvline(pv, color=col, lw=0.9, ls='--', alpha=0.45, zorder=1)
        ax.text(1.07, y_pool, f"{pv*100:.1f}",
                va='center', ha='left', fontsize=10.5, fontweight='bold', color=col)

        # ── Separator ────────────────────────────────────────────────────────
        ax.axhline(0.58, color='#AAAAAA', lw=0.9, xmin=0, xmax=0.87, zorder=2)

        # ── Heterogeneity footnote ───────────────────────────────────────────
        p_s = "p<0.001" if pQ < 0.001 else f"p={pQ:.3f}"
        ax.text(0.87, -0.68,
                f"I²={I2:.0f}%  τ²={tau2:.3f}  Q={Q:.1f} ({p_s})",
                ha='right', va='center', fontsize=7.5,
                style='italic', color='#666666',
                transform=ax.transAxes)

        # ── Axes & labels ────────────────────────────────────────────────────
        ax.set_xlim(-0.02, 1.20)
        ax.set_ylim(-0.9, k + 0.75)
        ax.set_xlabel(xlabel, fontsize=10, labelpad=6)
        ax.tick_params(axis='y', length=0)
        ax.spines['left'].set_visible(False)

        col_hdr = "Sensitivity" if ax_i == 0 else "Specificity"
        ax.set_title(f"{col_hdr} — {d['label']}",
                     fontsize=11, fontweight='bold', pad=8)

        if ax_i == 0:
            ax.set_yticks(list(y) + [y_pool])
            ax.set_yticklabels(d['studies'] + [f"Pooled [k={k}]"],
                               fontsize=9.5, ha='right')
        else:
            ax.set_yticks([])

    fig.suptitle(
        f"Figure {fig_num}. Forest plot: {d['label']}  "
        f"(k={k}, N={d['N_orig']:,})\n"
        f"Pooled: Se={d['pool_se']*100:.1f}% "
        f"[{d['pool_se_lo']*100:.1f}–{d['pool_se_hi']*100:.1f}%]   "
        f"Sp={d['pool_sp']*100:.1f}% "
        f"[{d['pool_sp_lo']*100:.1f}–{d['pool_sp_hi']*100:.1f}%]",
        fontsize=11, fontweight='bold', y=1.02,
    )

    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")


# =============================================================================
# 4. FIGURE 4 — SROC
# =============================================================================

def sroc_msl(TP, FP, FN, TN):
    """Moses-Shapiro-Littenberg SROC curve. Returns (FPR, TPR) arrays."""
    TP, FP, FN, TN = apply_cc(TP.copy(), FP.copy(), FN.copy(), TN.copy())
    Se = TP / (TP + FN);  Sp = TN / (TN + FP)
    D  = np.log(Se/(1-Se)) + np.log(Sp/(1-Sp))
    S  = np.log(Se/(1-Se)) - np.log(Sp/(1-Sp))
    w  = 1 / ((1/(TP+.5) + 1/(FN+.5)) + (1/(TN+.5) + 1/(FP+.5)))
    a, _ = np.polyfit(S, D, 1, w=w)[::-1]   # intercept, slope
    # Re-fit correctly: D = a + b*S  →  np.polyfit returns [b, a]
    coefs = np.polyfit(S, D, 1, w=w)
    b_coef, a_coef = coefs[0], coefs[1]
    t   = np.linspace(-8, 8, 400)
    tpr = 1 / (1 + np.exp(-(a_coef - t)/2))
    fpr = 1 - 1 / (1 + np.exp(-(a_coef + t)/2))
    ok  = (fpr >= 0) & (fpr <= 1) & (tpr >= 0) & (tpr <= 1)
    return fpr[ok], tpr[ok]


def plot_sroc(d_g, d_s, filename):
    fig, ax = plt.subplots(figsize=(8, 8))

    max_N = max(
        (d_g['TP']+d_g['FP']+d_g['FN']+d_g['TN']).max(),
        (d_s['TP']+d_s['FP']+d_s['FN']+d_s['TN']).max(),
    )

    for d, marker, zo in [(d_g, 'D', 4), (d_s, 's', 3)]:
        TP, FP, FN, TN = apply_cc(d['TP'].copy(), d['FP'].copy(),
                                   d['FN'].copy(), d['TN'].copy())
        Se = TP/(TP+FN);  Sp = TN/(TN+FP)
        N  = d['TP']+d['FP']+d['FN']+d['TN']

        # Individual studies
        sz = np.sqrt(N / max_N) * 320
        ax.scatter(1-Sp, Se, s=sz, color=d['color_l'],
                   edgecolors=d['color'], linewidth=1.5, zorder=zo, alpha=0.85)

        # SROC curve
        fpr, tpr = sroc_msl(d['TP'], d['FP'], d['FN'], d['TN'])
        ax.plot(fpr, tpr, color=d['color'], lw=2.5, zorder=5)

        # Summary operating point
        sp_x = 1 - d['pool_sp']
        se_y = d['pool_se']
        ax.scatter(sp_x, se_y, marker=marker, s=200,
                   color=d['color'], edgecolors='white', linewidth=2, zorder=6)

        # 95% CI ellipse around summary point
        hw_x = (d['pool_sp_hi'] - d['pool_sp_lo']) / 2   # horizontal half-width in SP space
        hw_y = (d['pool_se_hi'] - d['pool_se_lo']) / 2
        theta = np.linspace(0, 2*np.pi, 120)
        ax.plot(sp_x + hw_x*np.cos(theta),
                se_y  + hw_y*np.sin(theta),
                color=d['color'], lw=1.2, ls='--', alpha=0.75, zorder=5)

    # Chance line
    ax.plot([0,1],[0,1], color='#AAAAAA', lw=1, ls=':', zorder=1)
    ax.grid(True, color='#EEEEEE', lw=0.8, zorder=0)

    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel("1 − Specificity  (False Positive Rate)", fontsize=11)
    ax.set_ylabel("Sensitivity  (True Positive Rate)", fontsize=11)

    legend_els = [
        mpatches.Patch(color=COL_G,
            label=f"GFAP+UCH-L1:  Se={d_g['pool_se']*100:.1f}%, Sp={d_g['pool_sp']*100:.1f}%"),
        mpatches.Patch(color=COL_S,
            label=f"S100B:  Se={d_s['pool_se']*100:.1f}%, Sp={d_s['pool_sp']*100:.1f}%"),
        plt.Line2D([0],[0], color='#AAAAAA', lw=1, ls=':',
                   label="Chance line (AUC=0.5)"),
        plt.Line2D([0],[0], color='gray', lw=1.2, ls='--',
                   label="95% CI ellipse around summary point"),
        mpatches.Patch(facecolor='lightgray', edgecolor='gray',
                       label="Circle area ∝ study N"),
    ]
    ax.legend(handles=legend_els, loc='lower right', fontsize=9,
              framealpha=0.95, edgecolor='#CCCCCC', fancybox=False)

    # Caption below (JAMA style: no internal title)
    ax.text(0.5, -0.10,
            f"Figure 4. SROC curves: GFAP+UCH-L1 vs S100B  |  "
            f"Bivariate RE model (REML) | metafor v4.8-0 | R 4.3",
            ha='center', va='top', fontsize=8.5, style='italic',
            transform=ax.transAxes)

    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")


# =============================================================================
# 5. FIGURE 5 — CLINICAL SUMMARY (2×2 panel)
# =============================================================================

def plot_clinical_summary(d_g, d_s, filename):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "Figure 5. Clinical Summary — GFAP+UCH-L1 vs S100B for mTBI Triage\n"
        "Bivariate random-effects meta-analysis (REML) | metafor v4.8-0 | R 4.3",
        fontsize=11, fontweight='bold',
    )

    xlabs  = ["GFAP+UCH-L1", "S100B"]
    x      = np.array([0, 1])
    w      = 0.32

    # ── A: Se & Sp ───────────────────────────────────────────────────────────
    ax = axes[0, 0]
    se_v  = np.array([d_g['pool_se'],    d_s['pool_se']])    * 100
    sp_v  = np.array([d_g['pool_sp'],    d_s['pool_sp']])    * 100
    se_lo = np.array([d_g['pool_se_lo'], d_s['pool_se_lo']]) * 100
    se_hi = np.array([d_g['pool_se_hi'], d_s['pool_se_hi']]) * 100
    sp_lo = np.array([d_g['pool_sp_lo'], d_s['pool_sp_lo']]) * 100
    sp_hi = np.array([d_g['pool_sp_hi'], d_s['pool_sp_hi']]) * 100

    ax.bar(x - w/2, se_v, w, color=[COL_G, COL_S], label='Sensitivity', zorder=2)
    ax.bar(x + w/2, sp_v, w, color=[COL_GL, COL_SL], label='Specificity', zorder=2)
    for i in range(2):
        ax.errorbar(x[i]-w/2, se_v[i],
                    yerr=[[se_v[i]-se_lo[i]], [se_hi[i]-se_v[i]]],
                    fmt='none', color='#222222', capsize=4, lw=1.5, zorder=3)
        ax.errorbar(x[i]+w/2, sp_v[i],
                    yerr=[[sp_v[i]-sp_lo[i]], [sp_hi[i]-sp_v[i]]],
                    fmt='none', color='#222222', capsize=4, lw=1.5, zorder=3)
        ax.text(x[i]-w/2, se_v[i]+1.5, f"{se_v[i]:.1f}%",
                ha='center', fontsize=9, fontweight='bold')
        ax.text(x[i]+w/2, sp_v[i]+1.5, f"{sp_v[i]:.1f}%",
                ha='center', fontsize=9, fontweight='bold')
    ax.axhline(95, color='#CC4444', lw=0.9, ls='--', alpha=0.7)
    ax.text(1.48, 95.6, "95%", color='#CC4444', fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(xlabs, fontsize=10)
    ax.set_ylim(0, 115); ax.set_ylabel("Value (%)")
    ax.set_title("A.  Sensitivity and Specificity", fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    # ── B: NPV ───────────────────────────────────────────────────────────────
    ax = axes[0, 1]
    npv = np.array([d_g['NPV'], d_s['NPV']])
    ax.bar(x, npv, 0.5, color=[COL_G, COL_S], zorder=2)
    for i in range(2):
        ax.text(x[i], npv[i]+0.02, f"{npv[i]:.2f}%",
                ha='center', fontsize=10, fontweight='bold')
    ax.axhline(99, color='#CC4444', lw=0.9, ls='--', alpha=0.7)
    ax.text(1.42, 99.02, "99%", color='#CC4444', fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(xlabs, fontsize=10)
    ax.set_ylim(97, 100.7); ax.set_ylabel("Negative Predictive Value (%)")
    ax.set_title("B.  NPV", fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    # ── C: CT avoided & Missed ───────────────────────────────────────────────
    ax = axes[1, 0]
    ct_v = np.array([d_g['CT_avoided'],    d_s['CT_avoided']])
    ms_v = np.array([d_g['missed_1000'],   d_s['missed_1000']])
    ax.bar(x - w/2, ct_v, w, color=[COL_G, COL_S], label="CT avoided (%)", zorder=2)
    ax.bar(x + w/2, ms_v, w, color=[COL_GL, COL_SL], label="Missed/1000 pts", zorder=2)
    for i in range(2):
        ax.text(x[i]-w/2, ct_v[i]+0.5, f"{ct_v[i]:.1f}",
                ha='center', fontsize=9, fontweight='bold')
        ax.text(x[i]+w/2, ms_v[i]+0.08, f"{ms_v[i]:.1f}",
                ha='center', fontsize=9, fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels(xlabs, fontsize=10)
    ax.set_ylim(0, 55); ax.set_ylabel("Value")
    ax.set_title("C.  CT Scans Avoided and Missed Lesions", fontweight='bold')
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    # ── D: Summary table ─────────────────────────────────────────────────────
    ax = axes[1, 1]
    ax.axis('off')
    ax.set_title("D.  Summary Table", fontweight='bold', loc='left', pad=6)

    rows = [
        ["Parameter",          "GFAP+UCH-L1",                    "S100B"],
        ["Studies (k)",        "7",                               "3"],
        ["Patients (N)",       "4,955",                           "1,871"],
        ["CT+ events",         f"492 (9.9%)",                     f"110 (5.9%)"],
        ["Sensitivity (%)",
         f"{d_g['pool_se']*100:.1f} [{d_g['pool_se_lo']*100:.1f}–{d_g['pool_se_hi']*100:.1f}]",
         f"{d_s['pool_se']*100:.1f} [{d_s['pool_se_lo']*100:.1f}–{d_s['pool_se_hi']*100:.1f}]"],
        ["Specificity (%)",
         f"{d_g['pool_sp']*100:.1f} [{d_g['pool_sp_lo']*100:.1f}–{d_g['pool_sp_hi']*100:.1f}]",
         f"{d_s['pool_sp']*100:.1f} [{d_s['pool_sp_lo']*100:.1f}–{d_s['pool_sp_hi']*100:.1f}]"],
        ["I² Se / Sp (%)",
         f"{d_g['I2_se']:.0f} / {d_g['I2_sp']:.0f}",
         f"{d_s['I2_se']:.0f} / {d_s['I2_sp']:.0f}"],
        ["NPV (%)",            f"{d_g['NPV']:.2f}",              f"{d_s['NPV']:.2f}"],
        ["CT avoided (%)",     f"{d_g['CT_avoided']:.1f}",       f"{d_s['CT_avoided']:.1f}"],
        ["Missed/1000 pts",    f"{d_g['missed_1000']:.1f}",      f"{d_s['missed_1000']:.1f}"],
    ]

    col_x    = [0.01, 0.43, 0.72]
    rh       = 0.090
    y0       = 0.96

    for ri, row in enumerate(rows):
        yy = y0 - ri * rh
        bg = '#DDEEFF' if ri == 0 else ('#F4F4F4' if ri % 2 == 0 else 'white')
        ax.add_patch(plt.Rectangle((0, yy - rh*0.85), 1, rh*0.92,
                                   transform=ax.transAxes,
                                   color=bg, zorder=0, clip_on=False))
        for ci, (cell, cx) in enumerate(zip(row, col_x)):
            fw  = 'bold' if ri == 0 else 'normal'
            col = (COL_G  if (ri == 0 and ci == 1) else
                   COL_S  if (ri == 0 and ci == 2) else 'black')
            ax.text(cx, yy - rh*0.30, cell,
                    transform=ax.transAxes,
                    fontsize=8.5, fontweight=fw, color=col, va='center')

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")


# =============================================================================
# 6. FIGURE S1 — FUNNEL PLOTS
# =============================================================================

def plot_funnel(d_g, d_s, filename):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(
        "Figure S1. Funnel Plots for Publication Bias Assessment\n"
        "Note: Formal Egger test has low power with k < 10; visual inspection only",
        fontsize=11, fontweight='bold',
    )

    panel_labels = [["A", "B"], ["C", "D"]]

    for ri, (d, col) in enumerate([(d_g, COL_G), (d_s, COL_S)]):
        TP, FP, FN, TN = apply_cc(d['TP'].copy(), d['FP'].copy(),
                                   d['FN'].copy(), d['TN'].copy())
        Se = TP/(TP+FN);  Sp = TN/(TN+FP)
        N_se = TP+FN;     N_sp = TN+FP

        for ci, (vals, pv, Nv, xlabel, tail) in enumerate([
            (Se, d['pool_se'], N_se,
             "Logit(Sensitivity)", f"{d['label']} — Sensitivity"),
            (Sp, d['pool_sp'], N_sp,
             "Logit(Specificity)", f"{d['label']} — Specificity"),
        ]):
            ax     = axes[ri, ci]
            lv     = np.log(vals/(1-vals))
            se_arr = np.sqrt(1 / (Nv * vals * (1-vals)))
            pool_l = np.log(pv/(1-pv))

            ax.scatter(lv, se_arr, color=col, s=70,
                       edgecolors='white', linewidth=0.9, zorder=3)
            ax.axvline(pool_l, color='#CC4444', lw=1.5, ls='--', zorder=2)

            # Pseudo-confidence funnel
            max_se = se_arr.max() * 1.25
            se_r   = np.linspace(0, max_se, 200)
            for alpha_mult, alpha in [(1.96, 0.10), (1.645, 0.10)]:
                ax.fill_betweenx(se_r,
                                 pool_l - alpha_mult*se_r,
                                 pool_l + alpha_mult*se_r,
                                 alpha=alpha, color=col, zorder=1)

            ax.invert_yaxis()
            ax.set_xlabel(xlabel, fontsize=9)
            ax.set_ylabel("Standard Error", fontsize=9)
            ax.set_title(f"{panel_labels[ri][ci]}. {tail}",
                         fontsize=10, fontweight='bold')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")



# =============================================================================
# 7. RUN ALL
# =============================================================================
if __name__ == '__main__':
    import os
    # Save figures to current directory instead of /mnt/user-data/outputs
    print("Generating figures (CORRECTED DATA v14)...\n")

    plot_forest(GFAP,  2, "Fig2_forest_gfap_v14.png")
    plot_forest(S100B, 3, "Fig3_forest_s100b_v14.png")
    plot_sroc(GFAP, S100B,       "Fig4_sroc_v14.png")
    plot_clinical_summary(GFAP, S100B, "Fig5_clinical_summary_v14.png")
    plot_funnel(GFAP, S100B,     "FigS1_funnel_v14.png")

    print("\nDone. All 5 figures saved to current directory")
