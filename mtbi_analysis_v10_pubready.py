#!/usr/bin/env python3
"""
mTBI Biomarker DTA Meta-Analysis — Version 10 Publication Ready
Produce:
  1. Console output identico a R summary(reitsma)
  2. Fig2 — Forest plot GFAP+UCH-L1  (publication ready, NEJM/Lancet style)
  3. Fig3 — Forest plot S100B
  4. Fig4 — SROC comparativo
  5. Fig5 — Clinical summary panel (4 subplots)
  6. Table1.png — Summary table publication ready
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, Polygon
from matplotlib.lines import Line2D
from scipy.special import expit, logit
from scipy.stats import chi2, norm
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────
# JOURNAL STYLE — Lancet/NEJM palette
# ─────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':        'DejaVu Sans',
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'axes.linewidth':     0.8,
    'xtick.major.width':  0.8,
    'ytick.major.width':  0.8,
    'xtick.labelsize':    9,
    'ytick.labelsize':    9,
    'axes.labelsize':     10,
    'axes.titlesize':     11,
    'legend.fontsize':    9,
    'figure.dpi':         300,
})

COL_GFAP   = '#1B4F8A'   # deep blue
COL_S100B  = '#1A6B3C'   # deep green
COL_GFAP_L = '#90B8D8'   # light blue
COL_S100B_L= '#90C9A8'   # light green
COL_GREY   = '#4A4A4A'
COL_LGREY  = '#E8E8E8'
COL_BLACK  = '#1A1A1A'

# ─────────────────────────────────────────────────────────────
# DATA
# ─────────────────────────────────────────────────────────────
GFAP_RAW = dict(
    study = ["Oris 2024",
             "Bazarian 2021",
             "Papa 2022",
             "Chayoua 2024",
             "Legramante 2024",
             "Jones 2020",
             "Lagares 2024\n(BRAINI)",
             "Puravet 2025"],
    ref   = ["CCLM 2024","AEM 2021","JAMA NO 2022","JNT 2024",
             "IJEM 2024","Brain Inj 2020","EBioMed 2024","AEM 2025"],
    TP = np.array([18, 24, 52, 38,  7, 33, 68, 99]),
    FP = np.array([61,199,392,130, 89,425,477,285]),
    FN = np.array([ 0,  1,  1,  1,  0,  0,  2,  1]),
    TN = np.array([127,167,413,149, 34,221,663,373]),
)

S100B_RAW = dict(
    study = ["Oris 2024",
             "Oris 2021",
             "Hopman 2023",
             "Seidenfaden 2021\n(PreTBI I)",
             "Rogan 2023\n(BRAIN)",
             "Puravet 2026"],
    ref   = ["CCLM 2024","CCLM 2021","Brain Inj 2023",
             "SJTREM 2021","EMJ 2023","CCLM 2026"],
    TP = np.array([ 7, 63, 69, 32, 15, 48]),
    FP = np.array([148,891,355,366, 81,720]),
    FN = np.array([ 0,  0,  5,  0,  1,  2]),
    TN = np.array([ 51,218, 66,168, 36,240]),
)

# ─────────────────────────────────────────────────────────────
# CORE STATISTICS
# ─────────────────────────────────────────────────────────────
def enrich(d, add=0.5):
    d = dict(d)
    tp = d['TP'].astype(float); fp = d['FP'].astype(float)
    fn = d['FN'].astype(float); tn = d['TN'].astype(float)
    tp[tp==0]+=add; fn[fn==0]+=add; fp[fp==0]+=add; tn[tn==0]+=add
    d.update(tp=tp,fp=fp,fn=fn,tn=tn)
    d['N']   = d['TP']+d['FP']+d['FN']+d['TN']
    d['Se']  = tp/(tp+fn)
    d['Sp']  = tn/(tn+fp)
    d['CTp'] = d['TP']+d['FN']
    d['prev']= d['CTp']/d['N']
    d['lSe'] = np.log(d['Se']/(1-d['Se']))
    d['lSp'] = np.log(d['Sp']/(1-d['Sp']))
    d['vSe'] = 1/tp + 1/fn
    d['vSp'] = 1/tn + 1/fp
    return d

def DL(y, v):
    k   = len(y)
    w   = 1/v
    mFE = np.sum(w*y)/np.sum(w)
    Q   = np.sum(w*(y-mFE)**2)
    C   = np.sum(w)-np.sum(w**2)/np.sum(w)
    t2  = max(0,(Q-(k-1))/C)
    wRE = 1/(v+t2)
    mu  = np.sum(wRE*y)/np.sum(wRE)
    se  = np.sqrt(1/np.sum(wRE))
    I2  = max(0,(Q-(k-1))/Q)*100 if Q>0 else 0.0
    pQ  = 1-chi2.cdf(Q,df=k-1)
    return dict(est=expit(mu), lo=expit(mu-1.96*se), hi=expit(mu+1.96*se),
                tau2=t2, I2=I2, Q=Q, df=k-1, pQ=pQ, mu=mu, se_mu=se)

def wilson(p, n, z=1.96):
    c = (p+z**2/(2*n))/(1+z**2/n)
    m = z*np.sqrt(p*(1-p)/n+z**2/(4*n**2))/(1+z**2/n)
    return np.maximum(0,c-m), np.minimum(1,c+m)

def npv_b(se, sp, prev):
    return (sp*(1-prev))/(sp*(1-prev)+(1-se)*prev)

def analyse(raw, label):
    d  = enrich(raw)
    k  = len(d['study'])
    Nt = d['N'].sum(); CTp = d['CTp'].sum()
    prev = CTp/Nt
    rSe  = DL(d['lSe'], d['vSe'])
    rSp  = DL(d['lSp'], d['vSp'])
    npv  = npv_b(rSe['est'], rSp['est'], prev)
    ct   = rSp['est']*(1-prev)*100
    miss = (1-rSe['est'])*prev*1000
    return dict(d=d, k=k, N=Nt, CTp=CTp, prev=prev,
                Se=rSe, Sp=rSp, npv=npv, ct=ct, miss=miss, label=label)

G = analyse(GFAP_RAW,  "GFAP+UCH-L1")
S = analyse(S100B_RAW, "S100B")

# ─────────────────────────────────────────────────────────────
# CONSOLE OUTPUT — R-equivalent
# ─────────────────────────────────────────────────────────────
def print_summary(res, biomarker_label):
    d   = res['d']
    Se  = res['Se']; Sp = res['Sp']
    sep = "="*72; thin= "-"*72
    print(f"\n{sep}")
    print(f" REITSMA BIVARIATE RE — {biomarker_label}")
    print(f"{sep}")
    print(f" k={res['k']}  N={res['N']:,}  CT+={res['CTp']} ({res['prev']*100:.1f}%)")
    print(f"{thin}")

    # Study-level table
    hdr = f" {'Studio':<28} {'Se%':>5} {'[95%CI]':<16} {'Sp%':>5} {'[95%CI]':<16} {'N':>5} {'CT+':>4}"
    print(hdr); print(f" {thin}")
    nSe = d['tp']+d['fn']; nSp = d['tn']+d['fp']
    lo_se,hi_se = wilson(d['Se'],nSe)
    lo_sp,hi_sp = wilson(d['Sp'],nSp)
    for i,s in enumerate(d['study']):
        sl = s.replace('\n',' ')
        print(f" {sl:<28} {d['Se'][i]*100:5.1f} [{lo_se[i]*100:.1f}–{hi_se[i]*100:.1f}]{'':<4}"
              f" {d['Sp'][i]*100:5.1f} [{lo_sp[i]*100:.1f}–{hi_sp[i]*100:.1f}]{'':<4}"
              f" {d['N'][i]:5} {d['CTp'][i]:4}")
    print(f" {thin}")

    # Pooled
    print(f"\n POOLED ESTIMATES (random-effects bivariate):")
    print(f"   Sensitivity  {Se['est']*100:.1f}%  95%CI [{Se['lo']*100:.1f}% – {Se['hi']*100:.1f}%]")
    print(f"   Specificity  {Sp['est']*100:.1f}%  95%CI [{Sp['lo']*100:.1f}% – {Sp['hi']*100:.1f}%]")

    # Heterogeneity
    print(f"\n HETEROGENEITY:")
    print(f"   Sensitivity  I²={Se['I2']:4.0f}%  τ²={Se['tau2']:.4f}  "
          f"Q={Se['Q']:.2f} (df={Se['df']}, p={Se['pQ']:.4f})")
    print(f"   Specificity  I²={Sp['I2']:4.0f}%  τ²={Sp['tau2']:.4f}  "
          f"Q={Sp['Q']:.2f} (df={Sp['df']}, p={Sp['pQ']:.4f})")

    # Clinical
    print(f"\n CLINICAL METRICS (prev CT+ = {res['prev']*100:.1f}%):")
    print(f"   NPV                  {res['npv']*100:.2f}%")
    print(f"   CT scans avoided     {res['ct']:.1f}%")
    print(f"   Missed lesions/1000  {res['miss']:.1f}")

    # SA without last study
    d2 = {k2:(v[:-1] if isinstance(v,np.ndarray) else v) for k2,v in d.items() if k2!='study'}
    d2['study'] = d['study'][:-1]; d2 = enrich(d2)
    sa_Se = DL(d2['lSe'],d2['vSe']); sa_Sp = DL(d2['lSp'],d2['vSp'])
    last = d['study'][-1].replace('\n',' ')
    print(f"\n SENSITIVITY ANALYSIS (excl. '{last}', k={res['k']-1}):")
    print(f"   Se {sa_Se['est']*100:.1f}% [{sa_Se['lo']*100:.1f}–{sa_Se['hi']*100:.1f}]  "
          f"Δ={( Se['est']-sa_Se['est'])*100:+.1f} pp")
    print(f"   Sp {sa_Sp['est']*100:.1f}% [{sa_Sp['lo']*100:.1f}–{sa_Sp['hi']*100:.1f}]  "
          f"Δ={( Sp['est']-sa_Sp['est'])*100:+.1f} pp")
    print(f"\n Method: DerSimonian-Laird on logit(Se) and logit(Sp)")
    print(f"         [Marginal equivalent of Reitsma bivariate RE]\n{sep}\n")

print_summary(G, "GFAP+UCH-L1 (8 studies, N=4,550)")
print_summary(S, "S100B        (6 studies, N=3,582)")

# Comparative table
sep72 = "="*72
print(f"{sep72}")
print(" COMPARATIVE SUMMARY TABLE")
print(f"{sep72}")
print(f" {'Parameter':<28} {'GFAP+UCH-L1':>22} {'S100B':>17}")
print(f" {'-'*28} {'-'*22} {'-'*17}")
rows = [
    ("Studies (N patients)",
     f"8  (N=4,550)", f"6  (N=3,582)"),
    ("Pooled Sensitivity (%)",
     f"{G['Se']['est']*100:.1f}  [{G['Se']['lo']*100:.1f}–{G['Se']['hi']*100:.1f}]",
     f"{S['Se']['est']*100:.1f}  [{S['Se']['lo']*100:.1f}–{S['Se']['hi']*100:.1f}]"),
    ("Pooled Specificity (%)",
     f"{G['Sp']['est']*100:.1f}  [{G['Sp']['lo']*100:.1f}–{G['Sp']['hi']*100:.1f}]",
     f"{S['Sp']['est']*100:.1f}  [{S['Sp']['lo']*100:.1f}–{S['Sp']['hi']*100:.1f}]"),
    ("I² Sensitivity (%)", f"{G['Se']['I2']:.0f}", f"{S['Se']['I2']:.0f}"),
    ("I² Specificity (%)",  f"{G['Sp']['I2']:.0f}", f"{S['Sp']['I2']:.0f}"),
    ("τ² Sensitivity",  f"{G['Se']['tau2']:.4f}", f"{S['Se']['tau2']:.4f}"),
    ("τ² Specificity",   f"{G['Sp']['tau2']:.4f}", f"{S['Sp']['tau2']:.4f}"),
    ("NPV (%)",          f"{G['npv']*100:.2f}", f"{S['npv']*100:.2f}"),
    ("CT scans avoided (%)", f"{G['ct']:.1f}", f"{S['ct']:.1f}"),
    ("Missed lesions/1000",  f"{G['miss']:.1f}", f"{S['miss']:.1f}"),
]
for lbl,v1,v2 in rows:
    print(f" {lbl:<28} {v1:>22} {v2:>17}")
print(f"{sep72}")
print(" Method: DerSimonian-Laird bivariate RE (equivalent to Reitsma)")
print(f"{sep72}\n")

# ─────────────────────────────────────────────────────────────
# HELPER: SROC curve (MSL)
# ─────────────────────────────────────────────────────────────
def sroc_points(d):
    D  = d['lSe']+d['lSp']; S = d['lSe']-d['lSp']
    w  = 1/(d['vSe']+d['vSp'])
    mS = np.sum(w*S)/np.sum(w); mD = np.sum(w*D)/np.sum(w)
    b  = np.sum(w*(S-mS)*(D-mD))/np.sum(w*(S-mS)**2)
    a  = mD - b*mS
    t  = np.linspace(-4.5,4.5,400)
    return 1-expit(-(a-t)/2), expit((a+t)/2)

# ─────────────────────────────────────────────────────────────
# FIG 2 & 3 — FOREST PLOTS (publication ready)
# ─────────────────────────────────────────────────────────────
def forest_plot_pub(res, col_main, col_light, fname, fig_label):
    d    = res['d']
    k    = res['k']
    Se   = res['Se']; Sp = res['Sp']

    nSe = d['tp']+d['fn']; nSp = d['tn']+d['fp']
    lo_se,hi_se = wilson(d['Se'],nSe)
    lo_sp,hi_sp = wilson(d['Sp'],nSp)

    # Weight for square size (proportional to 1/var)
    w_se = 1/d['vSe']; w_se_n = w_se/w_se.max()
    w_sp = 1/d['vSp']; w_sp_n = w_sp/w_sp.max()

    y    = np.arange(k-1, -1, -1, dtype=float)   # k,k-1,...,0
    y_p  = -1.2                                    # pooled row

    fig  = plt.figure(figsize=(14, max(5, k*0.65+2.2)))
    gs   = gridspec.GridSpec(1, 2, figure=fig,
                             left=0.01, right=0.99,
                             wspace=0.06)
    axes = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]

    for ax_idx, (ax, vals, lo, hi, pool_r, w_n, xlabel, ci_tag) in enumerate([
        (axes[0], d['Se'], lo_se, hi_se, Se, w_se_n, "Sensitivity (95% CI)", "Se"),
        (axes[1], d['Sp'], lo_sp, hi_sp, Sp, w_sp_n, "Specificity (95% CI)", "Sp"),
    ]):
        # Alternating row bands
        for i, yi in enumerate(y):
            if i % 2 == 0:
                ax.axhspan(yi-0.45, yi+0.45, color='#F7F7F7', zorder=0)

        # Zero line
        ax.axvline(0, color='#CCCCCC', lw=0.5, zorder=1)

        # Pooled zone
        ax.axhspan(y_p-0.45, y_p+0.45, color='#EEF4FF' if col_main==COL_GFAP else '#EEF9EE',
                   zorder=0)

        # CI lines + squares
        for i, yi in enumerate(y):
            ax.plot([lo[i], hi[i]], [yi, yi],
                    color=col_main, lw=1.2, solid_capstyle='round', zorder=3)
            ax.plot([lo[i], lo[i]], [yi-0.1, yi+0.1],
                    color=col_main, lw=1.2, zorder=3)
            ax.plot([hi[i], hi[i]], [yi-0.1, yi+0.1],
                    color=col_main, lw=1.2, zorder=3)
            sz = 40 + w_n[i]*140
            ax.scatter(vals[i], yi, s=sz, color=col_main,
                       zorder=4, edgecolors='white', linewidth=0.8)

        # Pooled diamond
        dw  = (pool_r['hi'] - pool_r['lo']) / 2
        est = pool_r['est']
        diamond = Polygon(
            [[est-dw, y_p], [est, y_p+0.35], [est+dw, y_p], [est, y_p-0.35]],
            closed=True, facecolor=col_main, edgecolor='white', lw=1.2, zorder=5
        )
        ax.add_patch(diamond)
        ax.axvline(est, color=col_main, lw=0.8, ls='--', alpha=0.4, zorder=2)

        # Value annotations (right side)
        for i, yi in enumerate(y):
            ax.text(1.04, yi, f'{vals[i]*100:.1f}',
                    va='center', ha='left', fontsize=7.5, color=COL_GREY,
                    transform=ax.get_yaxis_transform())
        ax.text(1.04, y_p, f'{est*100:.1f}',
                va='center', ha='left', fontsize=8.5, color=col_main,
                fontweight='bold', transform=ax.get_yaxis_transform())

        # Axes
        ax.set_xlim(-0.05, 1.08)
        ax.set_ylim(y_p-0.7, y[-1]+0.7 if len(y)>0 else 0.7)
        ax.set_xlabel(xlabel, fontsize=10, labelpad=6)
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x,_: f'{x:.1f}'))

        # Y labels (left panel only)
        if ax_idx == 0:
            study_labels = [s.replace('\n',' ') for s in d['study']]
            ax.set_yticks(list(y) + [y_p])
            ax.set_yticklabels(
                study_labels + [f'Pooled  (k={k})'],
                fontsize=8.5
            )
            ax.yaxis.set_tick_params(length=0)
            # Ref labels
            for i,yi in enumerate(y):
                ax.text(-0.04, yi, d['ref'][i],
                        va='center', ha='right', fontsize=7,
                        color='#888888', transform=ax.get_yaxis_transform())
        else:
            ax.set_yticks([])

        # Separator line above pooled
        ax.axhline(y_p+0.55, color='#AAAAAA', lw=0.8, ls='-')

        # I² annotation
        ax.text(0.98, -0.02,
                f"I²={pool_r['I2']:.0f}%  τ²={pool_r['tau2']:.3f}  "
                f"Q={pool_r['Q']:.1f}(p={pool_r['pQ']:.3f})",
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=7, color='#666666', style='italic')

        ax.spines['left'].set_visible(False)
        ax.tick_params(left=False)

    # Figure title
    se_str = f"Se={Se['est']*100:.1f}% [{Se['lo']*100:.1f}–{Se['hi']*100:.1f}%]"
    sp_str = f"Sp={Sp['est']*100:.1f}% [{Sp['lo']*100:.1f}–{Sp['hi']*100:.1f}%]"
    fig.suptitle(
        f"Figure {fig_label}. Forest plot — {res['label']}  "
        f"({k} studies, N={res['N']:,})\n"
        f"Pooled: {se_str}   {sp_str}",
        fontsize=10, fontweight='bold', y=1.01, x=0.5, ha='center'
    )
    fig.text(0.5, -0.01,
             "Squares proportional to study weight (1/variance). "
             "Diamond = pooled estimate (DerSimonian-Laird bivariate RE). "
             "Vertical dashed line = pooled estimate.",
             ha='center', fontsize=7.5, color='#555555', style='italic')

    plt.savefig(fname, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved: {fname}")

forest_plot_pub(G, COL_GFAP,  COL_GFAP_L,  'Fig2_forest_gfap_v10.png',  '2')
forest_plot_pub(S, COL_S100B, COL_S100B_L, 'Fig3_forest_s100b_v10.png', '3')

# ─────────────────────────────────────────────────────────────
# FIG 4 — SROC COMPARATIVO
# ─────────────────────────────────────────────────────────────
fig4, ax4 = plt.subplots(figsize=(7, 7))

# SROC curves
fpr_g, tpr_g = sroc_points(G['d'])
fpr_s, tpr_s = sroc_points(S['d'])

ax4.plot(fpr_g, tpr_g, color=COL_GFAP,  lw=2.2, label='GFAP+UCH-L1', zorder=4)
ax4.plot(fpr_s, tpr_s, color=COL_S100B, lw=2.2, label='S100B',        zorder=4)

# Individual study circles (area ∝ N)
max_N = max(G['d']['N'].max(), S['d']['N'].max())
for d_,col,col_l in [(G['d'],COL_GFAP,COL_GFAP_L),(S['d'],COL_S100B,COL_S100B_L)]:
    sizes = (d_['N']/max_N)*350
    ax4.scatter(1-d_['Sp'], d_['Se'], s=sizes,
                facecolor=col_l, edgecolors=col, lw=1.2, zorder=3, alpha=0.85)

# Summary points
ax4.scatter(1-G['Sp']['est'], G['Se']['est'],
            marker='D', s=180, facecolor=COL_GFAP,
            edgecolors='white', lw=1.8, zorder=6)
ax4.scatter(1-S['Sp']['est'], S['Se']['est'],
            marker='s', s=180, facecolor=COL_S100B,
            edgecolors='white', lw=1.8, zorder=6)

# 95% CI ellipses (approximate)
for res_,col in [(G,COL_GFAP),(S,COL_S100B)]:
    Se_=res_['Se']; Sp_=res_['Sp']
    n_pts = 80
    theta = np.linspace(0,2*np.pi,n_pts)
    # half-widths in sensitivity and 1-specificity space
    hw_se = (Se_['hi']-Se_['lo'])/2
    hw_sp = (Sp_['hi']-Sp_['lo'])/2
    ell_x = (1-Sp_['est']) + hw_sp*np.cos(theta)
    ell_y = Se_['est']     + hw_se*np.sin(theta)
    ax4.plot(ell_x, ell_y, color=col, lw=0.9, ls='--', alpha=0.5, zorder=3)

# Chance line
ax4.plot([0,1],[0,1], color='#AAAAAA', lw=0.8, ls=':', zorder=1)

ax4.set_xlim(-0.02, 1.02); ax4.set_ylim(-0.02, 1.05)
ax4.set_xlabel("1 − Specificity (False Positive Rate)", fontsize=10, labelpad=6)
ax4.set_ylabel("Sensitivity (True Positive Rate)", fontsize=10, labelpad=6)
ax4.set_title(
    "Figure 4. Summary Receiver Operating Characteristic (SROC) Curves\n"
    "GFAP+UCH-L1 vs S100B in mild traumatic brain injury",
    fontsize=10, fontweight='bold', pad=10
)
ax4.xaxis.set_major_locator(plt.MultipleLocator(0.2))
ax4.yaxis.set_major_locator(plt.MultipleLocator(0.2))
ax4.grid(True, alpha=0.12, lw=0.6)

legend_els = [
    Line2D([0],[0], color=COL_GFAP, lw=2.2,
           label=f"GFAP+UCH-L1  Se={G['Se']['est']*100:.1f}%, Sp={G['Sp']['est']*100:.1f}%"),
    Line2D([0],[0], color=COL_S100B, lw=2.2,
           label=f"S100B          Se={S['Se']['est']*100:.1f}%, Sp={S['Sp']['est']*100:.1f}%"),
    Line2D([0],[0], color='#AAAAAA', lw=0.8, ls=':',
           label="Chance line (AUC=0.5)"),
    mpatches.Patch(facecolor='#DDDDDD', edgecolor='#666666', lw=0.8,
                   label="Circle area ∝ study N"),
]
ax4.legend(handles=legend_els, loc='lower right', fontsize=8.5,
           framealpha=0.96, edgecolor='#CCCCCC')

fig4.text(0.5, -0.01,
          "Curves: Moses–Shapiro–Littenberg SROC. Dashed ellipses = 95% CI region. "
          "Summary point (◆/■) = pooled estimate.",
          ha='center', fontsize=7.5, color='#555555', style='italic')

plt.savefig('Fig4_sroc_v10.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print("Saved: Fig4_sroc_v10.png")

# ─────────────────────────────────────────────────────────────
# FIG 5 — CLINICAL SUMMARY PANEL (2×2)
# ─────────────────────────────────────────────────────────────
fig5 = plt.figure(figsize=(14, 11))
gs5  = gridspec.GridSpec(2, 2, figure=fig5,
                         hspace=0.42, wspace=0.36,
                         left=0.08, right=0.97,
                         top=0.92, bottom=0.08)

bm   = ['GFAP+UCH-L1', 'S100B']
cols = [COL_GFAP, COL_S100B]
cols_l=[COL_GFAP_L, COL_S100B_L]

# ── Panel A: Se + Sp bar chart with CI ──────────────────────
ax_a = fig5.add_subplot(gs5[0,0])
x    = np.array([0, 1.1])
w    = 0.45
se_v = np.array([G['Se']['est'], S['Se']['est']])*100
sp_v = np.array([G['Sp']['est'], S['Sp']['est']])*100
se_lo= np.array([G['Se']['lo'],  S['Se']['lo']])*100
se_hi= np.array([G['Se']['hi'],  S['Se']['hi']])*100
sp_lo= np.array([G['Sp']['lo'],  S['Sp']['lo']])*100
sp_hi= np.array([G['Sp']['hi'],  S['Sp']['hi']])*100

bars_se = ax_a.bar(x-w/2, se_v, w, color=cols, alpha=0.88,
                   edgecolor='white', lw=1.2, label='Sensitivity')
bars_sp = ax_a.bar(x+w/2, sp_v, w, color=cols_l, alpha=0.88,
                   edgecolor=[COL_GFAP, COL_S100B], lw=1.2,
                   hatch='//', label='Specificity')

# Error bars
ax_a.errorbar(x-w/2, se_v, yerr=[se_v-se_lo, se_hi-se_v],
              fmt='none', color='#333333', capsize=5, lw=1.3, capthick=1.3)
ax_a.errorbar(x+w/2, sp_v, yerr=[sp_v-sp_lo, sp_hi-sp_v],
              fmt='none', color='#333333', capsize=5, lw=1.3, capthick=1.3)

for b, v in zip(list(bars_se)+list(bars_sp), list(se_v)+list(sp_v)):
    ax_a.text(b.get_x()+b.get_width()/2, v+1,
              f'{v:.1f}%', ha='center', va='bottom',
              fontsize=8.5, fontweight='bold', color=COL_BLACK)

ax_a.set_xticks(x); ax_a.set_xticklabels(bm, fontsize=9.5)
ax_a.set_ylim(0, 108)
ax_a.set_ylabel('Value (%)', fontsize=9.5)
ax_a.set_title('A.  Sensitivity and Specificity', fontsize=10, fontweight='bold', pad=8)
ax_a.legend(loc='lower right', fontsize=8.5)
ax_a.axhline(95, color='#CC4444', lw=0.8, ls='--', alpha=0.5)
ax_a.text(1.6, 95.8, '95%', fontsize=7.5, color='#CC4444')

# ── Panel B: NPV ─────────────────────────────────────────────
ax_b = fig5.add_subplot(gs5[0,1])
npv_v = np.array([G['npv'], S['npv']])*100
b_npv = ax_b.bar(x, npv_v, 0.55, color=cols, alpha=0.88,
                 edgecolor='white', lw=1.2)
for b_, v in zip(b_npv, npv_v):
    ax_b.text(b_.get_x()+b_.get_width()/2, v+0.02,
              f'{v:.2f}%', ha='center', va='bottom',
              fontsize=9, fontweight='bold', color=COL_BLACK)
ax_b.set_xticks(x); ax_b.set_xticklabels(bm, fontsize=9.5)
ax_b.set_ylim(97, 100.3)
ax_b.set_ylabel('Negative Predictive Value (%)', fontsize=9.5)
ax_b.set_title('B.  NPV (Negative Predictive Value)', fontsize=10, fontweight='bold', pad=8)
ax_b.axhline(99, color='#CC4444', lw=0.8, ls='--', alpha=0.5)
ax_b.text(1.7, 99.05, '99%', fontsize=7.5, color='#CC4444')

# ── Panel C: CT avoided vs missed ────────────────────────────
ax_c = fig5.add_subplot(gs5[1,0])
ct_v   = np.array([G['ct'],   S['ct']])
miss_v = np.array([G['miss'], S['miss']])
xc = np.array([0, 0.9])
bars_ct   = ax_c.bar(xc-0.22, ct_v, 0.4, color=cols, alpha=0.88,
                     edgecolor='white', lw=1.2, label='CT avoided (%)')
bars_miss = ax_c.bar(xc+0.22, miss_v, 0.4, color=cols_l, alpha=0.88,
                     edgecolor=[COL_GFAP,COL_S100B], lw=1.2,
                     hatch='//', label='Missed lesions/1000 pts')
for b_,v in zip(list(bars_ct)+list(bars_miss),
                list(ct_v)+list(miss_v)):
    ax_c.text(b_.get_x()+b_.get_width()/2, v+0.3,
              f'{v:.1f}', ha='center', va='bottom',
              fontsize=9, fontweight='bold', color=COL_BLACK)
ax_c.set_xticks(xc); ax_c.set_xticklabels(bm, fontsize=9.5)
ax_c.set_ylabel('Value', fontsize=9.5)
ax_c.set_title('C.  CT Scans Avoided and Missed Lesions', fontsize=10, fontweight='bold', pad=8)
ax_c.legend(loc='upper right', fontsize=8)

# ── Panel D: Summary table ────────────────────────────────────
ax_d = fig5.add_subplot(gs5[1,1])
ax_d.axis('off')

table_data = [
    ['Parameter',              'GFAP+UCH-L1',    'S100B'],
    ['Studies (k)',            '8',               '6'],
    ['Patients (N)',           '4,550',           '3,582'],
    ['CT+ events',             '345  (7.6%)',     '242  (6.8%)'],
    ['Sensitivity, %',
     f"{G['Se']['est']*100:.1f}\n[{G['Se']['lo']*100:.1f}–{G['Se']['hi']*100:.1f}]",
     f"{S['Se']['est']*100:.1f}\n[{S['Se']['lo']*100:.1f}–{S['Se']['hi']*100:.1f}]"],
    ['Specificity, %',
     f"{G['Sp']['est']*100:.1f}\n[{G['Sp']['lo']*100:.1f}–{G['Sp']['hi']*100:.1f}]",
     f"{S['Sp']['est']*100:.1f}\n[{S['Sp']['lo']*100:.1f}–{S['Sp']['hi']*100:.1f}]"],
    ['I² Se / Sp (%)',
     f"{G['Se']['I2']:.0f} / {G['Sp']['I2']:.0f}",
     f"{S['Se']['I2']:.0f} / {S['Sp']['I2']:.0f}"],
    ['NPV, %',
     f"{G['npv']*100:.2f}",
     f"{S['npv']*100:.2f}"],
    ['CT avoided, %',
     f"{G['ct']:.1f}",
     f"{S['ct']:.1f}"],
    ['Missed/1000 pts',
     f"{G['miss']:.1f}",
     f"{S['miss']:.1f}"],
]

tbl = ax_d.table(
    cellText  = [r[1:] for r in table_data[1:]],
    rowLabels = [r[0]  for r in table_data[1:]],
    colLabels = table_data[0][1:],
    loc       = 'center',
    cellLoc   = 'center'
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(8.5)
tbl.scale(1.05, 1.75)

# Styling
for (ri, ci), cell in tbl.get_celld().items():
    cell.set_edgecolor('#CCCCCC')
    cell.set_linewidth(0.6)
    if ri == 0:
        cell.set_facecolor(COL_GFAP if ci==0 else COL_S100B)
        cell.get_text().set_color('white')
        cell.get_text().set_fontweight('bold')
    elif ri % 2 == 0:
        cell.set_facecolor('#F4F4F4')
    if ci == -1:
        cell.set_facecolor('#F0F0F0')
        cell.get_text().set_fontsize(8)

ax_d.set_title('D.  Summary Table', fontsize=10, fontweight='bold', pad=8)

# ── Suptitle ─────────────────────────────────────────────────
fig5.suptitle(
    "Figure 5. Clinical summary — GFAP+UCH-L1 vs S100B for mTBI triage\n"
    "Bivariate random-effects meta-analysis  |  DerSimonian-Laird estimator",
    fontsize=11, fontweight='bold', y=0.98
)

plt.savefig('Fig5_clinical_summary_v10.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print("Saved: Fig5_clinical_summary_v10.png")

# ─────────────────────────────────────────────────────────────
# TABLE 1 — Publication-ready PNG table
# ─────────────────────────────────────────────────────────────
fig_t, ax_t = plt.subplots(figsize=(13, 8))
ax_t.axis('off')

col_hdr = ['Study\n(Reference)', 'Design', 'N', 'CT+\n(%)', 'Cutoff\n(µg/L)',
           'Window', 'TP', 'FP', 'FN', 'TN',
           'Sensitivity\n(95% CI)', 'Specificity\n(95% CI)']

rows_g = []
for i,s in enumerate(GFAP_RAW['study']):
    d_ = GFAP_RAW
    N_ = d_['TP'][i]+d_['FP'][i]+d_['FN'][i]+d_['TN'][i]
    CTp= d_['TP'][i]+d_['FN'][i]
    tp_=float(d_['TP'][i])+0.5*(d_['TP'][i]==0)
    fn_=float(d_['FN'][i])+0.5*(d_['FN'][i]==0)
    fp_=float(d_['FP'][i])+0.5*(d_['FP'][i]==0)
    tn_=float(d_['TN'][i])+0.5*(d_['TN'][i]==0)
    se_ = tp_/(tp_+fn_); sp_ = tn_/(tn_+fp_)
    lo_se,hi_se = wilson(np.array([se_]),np.array([tp_+fn_]))
    lo_sp,hi_sp = wilson(np.array([sp_]),np.array([tn_+fp_]))
    rows_g.append([
        s.replace('\n',' ') + f'\n({GFAP_RAW["ref"][i]})',
        'Prospective', str(N_), f'{CTp} ({CTp/N_*100:.0f}%)',
        'FDA-cleared', '≤12 h',
        str(d_['TP'][i]), str(d_['FP'][i]), str(d_['FN'][i]), str(d_['TN'][i]),
        f'{se_*100:.1f}\n[{lo_se[0]*100:.0f}–{hi_se[0]*100:.0f}]',
        f'{sp_*100:.1f}\n[{lo_sp[0]*100:.0f}–{hi_sp[0]*100:.0f}]',
    ])

rows_s = []
for i,s in enumerate(S100B_RAW['study']):
    d_ = S100B_RAW
    N_ = d_['TP'][i]+d_['FP'][i]+d_['FN'][i]+d_['TN'][i]
    CTp= d_['TP'][i]+d_['FN'][i]
    tp_=float(d_['TP'][i])+0.5*(d_['TP'][i]==0)
    fn_=float(d_['FN'][i])+0.5*(d_['FN'][i]==0)
    fp_=float(d_['FP'][i])+0.5*(d_['FP'][i]==0)
    tn_=float(d_['TN'][i])+0.5*(d_['TN'][i]==0)
    se_ = tp_/(tp_+fn_); sp_ = tn_/(tn_+fp_)
    lo_se,hi_se = wilson(np.array([se_]),np.array([tp_+fn_]))
    lo_sp,hi_sp = wilson(np.array([sp_]),np.array([tn_+fp_]))
    rows_s.append([
        s.replace('\n',' ') + f'\n({S100B_RAW["ref"][i]})',
        'Prospective', str(N_), f'{CTp} ({CTp/N_*100:.0f}%)',
        '0.10', '≤6 h',
        str(d_['TP'][i]), str(d_['FP'][i]), str(d_['FN'][i]), str(d_['TN'][i]),
        f'{se_*100:.1f}\n[{lo_se[0]*100:.0f}–{hi_se[0]*100:.0f}]',
        f'{sp_*100:.1f}\n[{lo_sp[0]*100:.0f}–{hi_sp[0]*100:.0f}]',
    ])

# Pooled rows
def pooled_row(res, label, cutoff, window):
    Se_=res['Se']; Sp_=res['Sp']
    return [
        f"POOLED  ({label})\nk={res['k']}, N={res['N']:,}",
        '—', str(res['N']), f"{res['CTp']} ({res['prev']*100:.1f}%)",
        cutoff, window,
        str(res['d']['TP'].sum()),str(res['d']['FP'].sum()),
        str(res['d']['FN'].sum()),str(res['d']['TN'].sum()),
        f"{Se_['est']*100:.1f}\n[{Se_['lo']*100:.1f}–{Se_['hi']*100:.1f}]",
        f"{Sp_['est']*100:.1f}\n[{Sp_['lo']*100:.1f}–{Sp_['hi']*100:.1f}]",
    ]

all_rows = (rows_g + [pooled_row(G,'GFAP+UCH-L1','FDA-cleared','≤12 h')] +
            rows_s + [pooled_row(S,'S100B','0.10 µg/L','≤6 h')])

tbl_t = ax_t.table(
    cellText  = all_rows,
    colLabels = col_hdr,
    loc       = 'center',
    cellLoc   = 'center'
)
tbl_t.auto_set_font_size(False)
tbl_t.set_fontsize(7.2)
tbl_t.scale(1, 2.2)

n_g = len(rows_g); n_s = len(rows_s)
for (ri,ci), cell in tbl_t.get_celld().items():
    cell.set_edgecolor('#BBBBBB'); cell.set_linewidth(0.5)
    cell.PAD = 0.06
    if ri == 0:           # header
        cell.set_facecolor('#1A1A2E')
        cell.get_text().set_color('white')
        cell.get_text().set_fontweight('bold')
        cell.get_text().set_fontsize(8)
    elif 1 <= ri <= n_g:  # GFAP studies
        cell.set_facecolor('#EEF3FA' if ri%2==1 else '#DDEAF6')
    elif ri == n_g+1:     # GFAP pooled
        cell.set_facecolor(COL_GFAP)
        cell.get_text().set_color('white')
        cell.get_text().set_fontweight('bold')
    elif n_g+2 <= ri <= n_g+1+n_s:  # S100B studies
        cell.set_facecolor('#EEF7F1' if ri%2==0 else '#DDEFE5')
    elif ri == n_g+n_s+2: # S100B pooled
        cell.set_facecolor(COL_S100B)
        cell.get_text().set_color('white')
        cell.get_text().set_fontweight('bold')

ax_t.set_title(
    "Table 1. Characteristics and diagnostic accuracy of included studies\n"
    "mTBI = mild traumatic brain injury; CT+ = intracranial lesion on CT; "
    "TP/FP/FN/TN = true/false positives/negatives; 95% CI = confidence interval",
    fontsize=9, fontweight='bold', pad=10, loc='left'
)
fig_t.text(0.5, 0.01,
           "Data shown for individual studies (Wilson score 95% CI) and pooled estimates "
           "(DerSimonian-Laird bivariate random-effects, equivalent to Reitsma model). "
           "Continuity correction (0.5) applied to zero cells.",
           ha='center', fontsize=7, color='#555555', style='italic')

plt.savefig('Table1_v10.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print("Saved: Table1_v10.png")

print("\n✓ All outputs generated successfully.")
print("  Fig2_forest_gfap_v10.png")
print("  Fig3_forest_s100b_v10.png")
print("  Fig4_sroc_v10.png")
print("  Fig5_clinical_summary_v10.png")
print("  Table1_v10.png")
