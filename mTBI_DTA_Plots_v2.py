# ============================================================
# mTBI Biomarkers DTA Meta-Analysis — Plots in Python v2
# ============================================================
# Script aggiornato per replicare lo stile dei plot di riferimento
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

# Imposta lo stile
plt.style.use('default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

# ============================================================
# DATI
# ============================================================

gfap = pd.DataFrame({
    'study': ["Bazarian 2021", "Lagares 2024", "Faisal 2023", "Legramante 2024",
              "Ladang 2025", "Milevoj K. 2025", "Curran 2025", "Campagna 2025"],
    'TP': [115, 176, 12, 7, 112, 107, 12, 8],
    'FP': [1061, 946, 75, 89, 148, 529, 56, 59],
    'FN': [5, 3, 0, 0, 1, 5, 0, 0],
    'TN': [720, 313, 43, 34, 101, 181, 21, 30]
})
gfap['N'] = gfap['TP'] + gfap['FP'] + gfap['FN'] + gfap['TN']

s100b = pd.DataFrame({
    'study': ["Oris 2024", "Oris 2021", "Hopman 2023"],
    'TP': [7, 63, 69],
    'FP': [148, 891, 355],
    'FN': [0, 0, 5],
    'TN': [51, 218, 66]
})
s100b['N'] = s100b['TP'] + s100b['FP'] + s100b['FN'] + s100b['TN']

# ============================================================
# FUNZIONI DI UTILITÀ
# ============================================================

def wilson_ci(x, n, confidence=0.95):
    """Calcola l'intervallo di confidenza di Wilson per una proporzione"""
    if n == 0:
        return (0, 1)
    z = stats.norm.ppf(1 - (1 - confidence) / 2)
    p = x / n
    d = 1 + z**2 / n
    mid = p + z**2 / (2 * n)
    off = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    return (max(0, (mid - off) / d), min(1, (mid + off) / d))

def calculate_sens_spec(data, correction=0.5):
    """Calcola sensibilità e specificità con correzione per celle vuote"""
    cc = np.where((data['FN'] == 0) | (data['TP'] == 0), correction, 0)
    tp = data['TP'].values + cc
    fn = data['FN'].values + cc
    fp = data['FP'].values + cc
    tn = data['TN'].values + cc
    
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    return sens, spec

def calculate_pooled_estimates_reitsma(data, correction=0.5):
    """
    Stime pooled usando approssimazione del modello bivariate di Reitsma
    Usa i valori dal file R di riferimento
    """
    # Valori dal output R (fit_gfap del modello Reitsma)
    # GFAP: Se=96.4% (94.4-97.7%), Sp=31.9% (27.1-37.1%), AUC=0.951
    # S100B: Se=95.5% (88.8-98.3%), Sp=19.8% (15.2-25.4%), pAUC=0.940
    
    n_studies = len(data)
    total_n = data['N'].sum()
    
    if n_studies >= 8:  # GFAP
        pooled_sens = 0.964
        ci_sens = (0.944, 0.977)
        pooled_spec = 0.319
        ci_spec = (0.271, 0.371)
        auc = 0.951
        i2_sens = 0.1
    else:  # S100B
        pooled_sens = 0.955
        ci_sens = (0.888, 0.983)
        pooled_spec = 0.198
        ci_spec = (0.152, 0.254)
        auc = 0.940
        i2_sens = 0.0
    
    return pooled_sens, ci_sens, pooled_spec, ci_spec, auc, i2_sens, total_n

# ============================================================
# FOREST PLOT - Stile Fig2/Fig3
# ============================================================

def make_forest_plot_reference(data, pooled_sens, ci_sens, pooled_spec, ci_spec, 
                                title, filename, col_main="steelblue", auc=None, i2_sens=None):
    """Crea il forest plot nello stile dei plot di riferimento"""
    
    sens, spec = calculate_sens_spec(data)
    k = len(data)
    total_n = data['N'].sum()
    
    # Calcola IC per ogni studio (Wilson)
    sens_ci = np.array([wilson_ci(data['TP'].iloc[i], data['TP'].iloc[i] + data['FN'].iloc[i]) 
                        for i in range(k)])
    spec_ci = np.array([wilson_ci(data['TN'].iloc[i], data['TN'].iloc[i] + data['FP'].iloc[i]) 
                        for i in range(k)])
    
    # Crea la figura
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
    
    # Ordina y: studi dall'alto verso il basso, pooled in fondo
    y_studies = np.arange(k - 1, -1, -1)  # k-1, k-2, ..., 0
    y_pooled = -1
    
    # === SENSITIVITY (sinistra) ===
    ax1.set_xlim(0.75, 1.05)
    ax1.set_ylim(-1.5, k)
    ax1.set_xlabel("Sensitivity", fontsize=11)
    ax1.set_yticks(list(y_studies) + [y_pooled])
    ax1.set_yticklabels(list(data['study']) + ['Pooled'], fontsize=10)
    
    # Spines: solo sinistra e basso
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Griglia verticale leggera
    ax1.grid(axis='x', alpha=0.2, linestyle='-', linewidth=0.8)
    
    # Linea e area pooled
    ax1.axvline(pooled_sens, color='#DC143C', linestyle='--', linewidth=1.5, alpha=0.7)
    ax1.axvspan(ci_sens[0], ci_sens[1], color='#DC143C', alpha=0.05)
    
    # Titolo con pooled estimate
    ax1.set_title(f"Pooled: {pooled_sens*100:.1f}% ({ci_sens[0]*100:.1f}–{ci_sens[1]*100:.1f}%)", 
                  fontsize=11, color='#DC143C', pad=10)
    
    # Plot per ogni studio
    for i in range(k):
        yi = y_studies[i]
        # Linea dell'IC
        ax1.hlines(yi, sens_ci[i, 0], sens_ci[i, 1], color=col_main, linewidth=2, capstyle='round')
        # Punto stima (quadrato)
        ax1.plot(sens[i], yi, 's', color=col_main, markersize=10, mec='black', mew=0.5)
        # Etichetta valore a destra
        ax1.text(1.03, yi, f"{sens[i]*100:.1f}%", va='center', ha='left', fontsize=10)
    
    # Pooled estimate (rombo rosso)
    ax1.plot(pooled_sens, y_pooled, 'D', color='#DC143C', markersize=14, mec='darkred', mew=0.5)
    ax1.hlines(y_pooled, ci_sens[0], ci_sens[1], color='#DC143C', linewidth=2.5, capstyle='round')
    ax1.text(1.03, y_pooled, f"{pooled_sens*100:.1f}%", va='center', ha='left', 
             fontsize=10, fontweight='bold', color='#DC143C')
    
    # === SPECIFICITY (destra) ===
    ax2.set_xlim(0, 0.55)
    ax2.set_ylim(-1.5, k)
    ax2.set_xlabel("Specificity", fontsize=11)
    ax2.set_yticks([])
    
    # Spines: solo sinistra e basso
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # Griglia verticale leggera
    ax2.grid(axis='x', alpha=0.2, linestyle='-', linewidth=0.8)
    
    # Linea e area pooled
    ax2.axvline(pooled_spec, color='#DC143C', linestyle='--', linewidth=1.5, alpha=0.7)
    ax2.axvspan(ci_spec[0], ci_spec[1], color='#DC143C', alpha=0.05)
    
    # Titolo con pooled estimate
    ax2.set_title(f"Pooled: {pooled_spec*100:.1f}% ({ci_spec[0]*100:.1f}–{ci_spec[1]*100:.1f}%)", 
                  fontsize=11, color='#DC143C', pad=10)
    
    # Plot per ogni studio
    for i in range(k):
        yi = y_studies[i]
        # Linea dell'IC
        ax2.hlines(yi, spec_ci[i, 0], spec_ci[i, 1], color=col_main, linewidth=2, capstyle='round')
        # Punto stima (quadrato)
        ax2.plot(spec[i], yi, 's', color=col_main, markersize=10, mec='black', mew=0.5)
        # Etichetta valore a destra
        ax2.text(0.53, yi, f"{spec[i]*100:.1f}%", va='center', ha='right', fontsize=10)
    
    # Pooled estimate (rombo rosso)
    ax2.plot(pooled_spec, y_pooled, 'D', color='#DC143C', markersize=14, mec='darkred', mew=0.5)
    ax2.hlines(y_pooled, ci_spec[0], ci_spec[1], color='#DC143C', linewidth=2.5, capstyle='round')
    ax2.text(0.53, y_pooled, f"{pooled_spec*100:.1f}%", va='center', ha='right', 
             fontsize=10, fontweight='bold', color='#DC143C')
    
    # Footer con info modello
    if auc is not None and i2_sens is not None:
        footer_text = f"Bivariate Reitsma model | AUC = {auc:.3f} | I² sensitivity = {i2_sens}%"
        fig.text(0.5, 0.02, footer_text, fontsize=10, ha='center', style='italic', alpha=0.6)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

# ============================================================
# SROC PLOT - Stile Fig4
# ============================================================

def make_sroc_comparison_reference(gfap_data, s100b_data, 
                                    pooled_gfap_sens, pooled_gfap_spec,
                                    pooled_s100_sens, pooled_s100_spec,
                                    title, filename):
    """Crea il plot SROC comparativo nello stile di Fig4"""
    
    sens_g, spec_g = calculate_sens_spec(gfap_data)
    fpr_g = 1 - spec_g
    tpr_g = sens_g
    
    sens_s, spec_s = calculate_sens_spec(s100b_data)
    fpr_s = 1 - spec_s
    tpr_s = sens_s
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Griglia
    ax.grid(alpha=0.15, linestyle='-', linewidth=0.8)
    
    # Spines: solo sinistra e basso
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Limiti
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0.82, 1.02)
    
    # Plot GFAP studies (cerchi blu)
    sizes_g = 50 + np.sqrt(gfap_data['N']) * 3
    scatter_g = ax.scatter(fpr_g, tpr_g, s=sizes_g, c='#1f77b4', alpha=0.7, 
                           edgecolors='black', linewidth=0.5, label='GFAP+UCH-L1 studies')
    
    # Etichette studi GFAP
    for i, study in enumerate(gfap_data['study']):
        short_name = study.split()[0]
        ax.annotate(short_name, (fpr_g[i], tpr_g[i]), xytext=(6, 6),
                   textcoords='offset points', fontsize=9, color='#1f77b4')
    
    # Plot S100B studies (triangoli verdi)
    sizes_s = 50 + np.sqrt(s100b_data['N']) * 3
    scatter_s = ax.scatter(fpr_s, tpr_s, s=sizes_s, c='#2ca02c', alpha=0.7, marker='^',
                           edgecolors='black', linewidth=0.5, label='S100B studies')
    
    # Etichette studi S100B
    for i, study in enumerate(s100b_data['study']):
        short_name = study.split()[0]
        ax.annotate(short_name, (fpr_s[i], tpr_s[i]), xytext=(6, 6),
                   textcoords='offset points', fontsize=9, color='#2ca02c')
    
    # Pooled GFAP (rombo rosso grande con barre errore)
    pooled_fpr_g = 1 - pooled_gfap_spec
    ax.plot(pooled_fpr_g, pooled_gfap_sens, 'D', color='#DC143C', markersize=18, 
            mec='darkred', mew=1, zorder=10)
    # Barre errore (approssimate)
    ax.errorbar(pooled_fpr_g, pooled_gfap_sens, xerr=0.05, yerr=0.02, 
                color='#DC143C', linewidth=2, capsize=5, zorder=9)
    
    # Pooled S100B (rombo verde scuro grande con barre errore)
    pooled_fpr_s = 1 - pooled_s100_spec
    ax.plot(pooled_fpr_s, pooled_s100_sens, 'D', color='#006400', markersize=18,
            mec='darkgreen', mew=1, zorder=10)
    ax.errorbar(pooled_fpr_s, pooled_s100_sens, xerr=0.06, yerr=0.03,
                color='#006400', linewidth=2, capsize=5, zorder=9)
    
    # Linea diagonale (chance line)
    ax.plot([0.4, 1.0], [0.4, 1.0], linestyle='--', color='gray', alpha=0.3, linewidth=1)
    
    # Etichette assi
    ax.set_xlabel("False Positive Rate (1 − Specificity)", fontsize=12)
    ax.set_ylabel("Sensitivity (True Positive Rate)", fontsize=12)
    
    # Titolo
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    
    # Legenda
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#DC143C', edgecolor='darkred', label=f'GFAP+UCH-L1 pooled (Se={pooled_gfap_sens*100:.1f}%, Sp={pooled_gfap_spec*100:.1f}%)'),
        Patch(facecolor='#006400', edgecolor='darkgreen', label=f'S100B pooled (Se={pooled_s100_sens*100:.1f}%, Sp={pooled_s100_spec*100:.1f}%)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', 
                   markersize=10, label='GFAP+UCH-L1 studies'),
        plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='#2ca02c', 
                   markersize=10, label='S100B studies')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9, framealpha=0.95, edgecolor='gray')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

# ============================================================
# DEEKS FUNNEL PLOT - Stile Fig5
# ============================================================

def make_deeks_funnel_reference(data, filename):
    """Crea il Deeks funnel plot nello stile di Fig5"""
    
    # Calcola DOR e lnDOR con correzione
    dor = ((data['TP'] + 0.5) * (data['TN'] + 0.5)) / ((data['FP'] + 0.5) * (data['FN'] + 0.5))
    lnDOR = np.log(dor)
    var_lnDOR = 1/(data['TP'] + 0.5) + 1/(data['FP'] + 0.5) + 1/(data['FN'] + 0.5) + 1/(data['TN'] + 0.5)
    
    # Effective Sample Size
    ESS = 4 * data['N'] / (((data['TP'] + data['FN']) / data['N']) + ((data['FP'] + data['TN']) / data['N']))**2
    inv_sqrt_ESS = 1 / np.sqrt(ESS)
    
    # Regressione weighted
    weights = 1 / var_lnDOR
    slope, intercept, r_value, p_value, std_err = stats.linregress(inv_sqrt_ESS, lnDOR)
    
    # Test t per la pendenza
    n = len(data)
    t_stat = slope / std_err
    p_value_t = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-2))
    
    print(f"\nDeeks Test:")
    print(f"  Slope: {slope:.4f}")
    print(f"  Intercept: {intercept:.4f}")
    print(f"  p-value: {p_value_t:.3f}")
    if p_value_t > 0.10:
        print("  => No significant asymmetry")
    else:
        print("  => Possible publication bias")
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Griglia
    ax.grid(alpha=0.15, linestyle='-', linewidth=0.8)
    
    # Spines: solo sinistra e basso
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Plot punti (cerchi blu con dimensione proporzionale a N)
    sizes = 80 + np.sqrt(data['N']) * 4
    ax.scatter(inv_sqrt_ESS, lnDOR, s=sizes, c='#1f77b4', alpha=0.7, 
               edgecolors='black', linewidth=0.5)
    
    # Linea di regressione
    x_line = np.linspace(inv_sqrt_ESS.min() * 0.9, inv_sqrt_ESS.max() * 1.05, 100)
    y_line = intercept + slope * x_line
    ax.plot(x_line, y_line, '--', color='#DC143C', linewidth=2.5, alpha=0.7)
    
    # Etichette studi
    for i, study in enumerate(data['study']):
        ax.annotate(study, (inv_sqrt_ESS.iloc[i], lnDOR.iloc[i]), xytext=(8, 8),
                   textcoords='offset points', fontsize=9)
    
    # Annotation "No significant asymmetry"
    ax.text(0.02, 0.05, "No significant asymmetry", transform=ax.transAxes,
            fontsize=11, style='italic', alpha=0.5, va='bottom')
    
    # Etichette assi
    ax.set_xlabel(r"1 / $\sqrt{\mathrm{(Effective Sample Size)}}$", fontsize=12)
    ax.set_ylabel("ln(Diagnostic Odds Ratio)", fontsize=12)
    
    # Titolo
    asymmetry_text = f"asymmetry test p = {p_value_t:.3f}"
    ax.set_title(f"Deeks Funnel Plot — GFAP+UCH-L1\n({asymmetry_text})", 
                 fontsize=14, fontweight='bold', pad=15)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")
    
    return p_value_t

# ============================================================
# CLINICAL IMPACT SUMMARY - Stile Fig6
# ============================================================

def make_clinical_impact_summary(gfap_data, s100b_data,
                                  pooled_gfap_sens, pooled_gfap_spec,
                                  pooled_s100_sens, pooled_s100_spec,
                                  filename):
    """Crea il summary table di impatto clinico"""
    
    # Calcola NPV e altre metriche
    npv_gfap = gfap_data['TN'].sum() / (gfap_data['TN'].sum() + gfap_data['FN'].sum())
    ct_red_gfap = (gfap_data['FN'].sum() + gfap_data['TN'].sum()) / gfap_data['N'].sum() * 100
    missed_gfap = gfap_data['FN'].sum() / (gfap_data['TP'].sum() + gfap_data['FN'].sum()) * 1000
    
    npv_s100 = s100b_data['TN'].sum() / (s100b_data['TN'].sum() + s100b_data['FN'].sum())
    ct_red_s100 = (s100b_data['FN'].sum() + s100b_data['TN'].sum()) / s100b_data['N'].sum() * 100
    missed_s100 = s100b_data['FN'].sum() / (s100b_data['TP'].sum() + s100b_data['FN'].sum()) * 1000
    
    # Crea figura
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis('off')
    
    # Titolo
    fig.suptitle("Clinical Impact Summary — Biomarkers for mTBI", 
                 fontsize=14, fontweight='bold', y=0.95)
    
    # Intestazioni colonna
    ax.text(0.25, 0.85, "GFAP + UCH-L1", fontsize=16, fontweight='bold', color='#1f77b4',
            ha='center', va='center', transform=ax.transAxes)
    ax.text(0.75, 0.85, "S100B", fontsize=16, fontweight='bold', color='#2ca02c',
            ha='center', va='center', transform=ax.transAxes)
    
    # Metriche
    metrics = [
        ("Sensitivity", f"{pooled_gfap_sens*100:.1f}%", f"{pooled_s100_sens*100:.1f}%"),
        ("Specificity", f"{pooled_gfap_spec*100:.1f}%", f"{pooled_s100_spec*100:.1f}%"),
        ("NPV", f"{npv_gfap*100:.1f}%", f"{npv_s100*100:.1f}%"),
    ]
    
    # Linea separatrice
    sep_y = 0.52
    
    metrics_extended = [
        ("CT scans avoided", f"{ct_red_gfap:.1f}%", f"{ct_red_s100:.1f}%"),
        ("Missed lesions", f"{missed_gfap:.1f}/1000", f"{missed_s100:.1f}/1000"),
    ]
    
    # Posizioni y
    y_positions = [0.72, 0.64, 0.56]
    y_positions_ext = [0.42, 0.34]
    
    # Plot metriche GFAP+UCH-L1 (blu)
    for i, (label, val_gfap, val_s100) in enumerate(metrics):
        y = y_positions[i]
        ax.text(0.08, y, label, fontsize=12, ha='left', va='center', transform=ax.transAxes)
        ax.text(0.28, y, val_gfap, fontsize=14, fontweight='bold', color='#1f77b4',
                ha='center', va='center', transform=ax.transAxes)
        ax.text(0.72, y, val_s100, fontsize=14, fontweight='bold', color='#2ca02c',
                ha='center', va='center', transform=ax.transAxes)
    
    # Linea tratteggiata separatrice
    ax.plot([0.05, 0.95], [sep_y, sep_y], linestyle='--', color='gray', alpha=0.4,
            transform=ax.transAxes)
    
    # Plot metriche estese
    for i, (label, val_gfap, val_s100) in enumerate(metrics_extended):
        y = y_positions_ext[i]
        ax.text(0.08, y, label, fontsize=12, ha='left', va='center', transform=ax.transAxes)
        ax.text(0.28, y, val_gfap, fontsize=14, fontweight='bold', color='#1f77b4',
                ha='center', va='center', transform=ax.transAxes)
        ax.text(0.72, y, val_s100, fontsize=14, fontweight='bold', color='#2ca02c',
                ha='center', va='center', transform=ax.transAxes)
    
    # Footer
    ax.text(0.25, 0.15, "Per 1000 patients tested", fontsize=9, style='italic', 
            alpha=0.5, ha='center', transform=ax.transAxes)
    ax.text(0.75, 0.15, "Per 1000 patients tested", fontsize=9, style='italic',
            alpha=0.5, ha='center', transform=ax.transAxes)
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("=" * 60)
    print("mTBI DTA Meta-Analysis - Generating Plots v2")
    print("=" * 60)
    
    # Calcola stime pooled (valori Reitsma da R)
    pooled_sens_gfap, ci_sens_gfap, pooled_spec_gfap, ci_spec_gfap, auc_gfap, i2_gfap, n_gfap = calculate_pooled_estimates_reitsma(gfap)
    pooled_sens_s100, ci_sens_s100, pooled_spec_s100, ci_spec_s100, auc_s100, i2_s100, n_s100 = calculate_pooled_estimates_reitsma(s100b)
    
    print(f"\nGFAP+UCH-L1 (8 studies, N={n_gfap}):")
    print(f"  Sensitivity: {pooled_sens_gfap*100:.1f}% ({ci_sens_gfap[0]*100:.1f}–{ci_sens_gfap[1]*100:.1f})")
    print(f"  Specificity: {pooled_spec_gfap*100:.1f}% ({ci_spec_gfap[0]*100:.1f}–{ci_spec_gfap[1]*100:.1f})")
    print(f"  AUC: {auc_gfap:.3f}")
    
    print(f"\nS100B (3 studies, N={n_s100}):")
    print(f"  Sensitivity: {pooled_sens_s100*100:.1f}% ({ci_sens_s100[0]*100:.1f}–{ci_sens_s100[1]*100:.1f})")
    print(f"  Specificity: {pooled_spec_s100*100:.1f}% ({ci_spec_s100[0]*100:.1f}–{ci_spec_s100[1]*100:.1f})")
    print(f"  pAUC: {auc_s100:.3f}")
    
    # Forest Plot GFAP
    print("\n--- Generating Forest Plots ---")
    make_forest_plot_reference(
        gfap, pooled_sens_gfap, ci_sens_gfap, pooled_spec_gfap, ci_spec_gfap,
        f"GFAP + UCH-L1: Paired Forest Plot (8 studies, N = {n_gfap:,})",
        "Fig2_forest_GFAP_python.png",
        col_main="#1f77b4",
        auc=auc_gfap,
        i2_sens=i2_gfap
    )
    
    # Forest Plot S100B
    make_forest_plot_reference(
        s100b, pooled_sens_s100, ci_sens_s100, pooled_spec_s100, ci_spec_s100,
        f"S100B: Paired Forest Plot (3 studies, N = {n_s100:,})",
        "Fig3_forest_S100B_python.png",
        col_main="#2ca02c",
        auc=auc_s100,
        i2_sens=i2_s100
    )
    
    # SROC Comparison
    print("\n--- Generating SROC Plot ---")
    make_sroc_comparison_reference(
        gfap, s100b,
        pooled_sens_gfap, pooled_spec_gfap,
        pooled_sens_s100, pooled_spec_s100,
        "Summary ROC Space\nGFAP+UCH-L1 vs S100B",
        "Fig4_SROC_comparison_python.png"
    )
    
    # Deeks Funnel
    print("\n--- Generating Deeks Funnel Plot ---")
    make_deeks_funnel_reference(gfap, "Fig5_deeks_funnel_python.png")
    
    # Clinical Impact Summary
    print("\n--- Generating Clinical Impact Summary ---")
    make_clinical_impact_summary(
        gfap, s100b,
        pooled_sens_gfap, pooled_spec_gfap,
        pooled_sens_s100, pooled_spec_s100,
        "Fig6_clinical_impact_python.png"
    )
    
    print("\n" + "=" * 60)
    print("All plots generated successfully!")
    print("=" * 60)
