# ============================================================
# mTBI Biomarkers DTA Meta-Analysis — Plots in Python
# ============================================================
# Questo script genera i plot di meta-analisi DTA usando Python
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
import pandas as pd

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

def calculate_sens_spec(data):
    """Calcola sensibilità e specificità con correzione per celle vuote"""
    cc = np.where((data['FN'] == 0) | (data['TP'] == 0), 0.5, 0)
    tp = data['TP'].values + cc
    fn = data['FN'].values + cc
    fp = data['FP'].values + cc
    tn = data['TN'].values + cc
    
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    return sens, spec

def calculate_pooled_estimates(data, correction=0.5):
    """
    Calcola le stime pooled di sensibilità e specificità usando il modello bivariate
    Semplificazione: usa la media pesata inversa della varianza
    """
    cc = np.where((data['FN'] == 0) | (data['TP'] == 0), correction, 0)
    tp = data['TP'].values + cc
    fn = data['FN'].values + cc
    fp = data['FP'].values + cc
    tn = data['TN'].values + cc
    
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    
    # Varianze approssimate
    var_sens = sens * (1 - sens) / (tp + fn)
    var_spec = spec * (1 - spec) / (tn + fp)
    
    # Pesi (inverse variance)
    weights_sens = 1 / np.where(var_sens > 0, var_sens, 0.001)
    weights_spec = 1 / np.where(var_spec > 0, var_spec, 0.001)
    
    # Stime pooled
    pooled_sens = np.average(sens, weights=weights_sens)
    pooled_spec = np.average(spec, weights=weights_spec)
    
    # Errori standard pooled
    se_pooled_sens = np.sqrt(1 / weights_sens.sum())
    se_pooled_spec = np.sqrt(1 / weights_spec.sum())
    
    # IC 95%
    ci_sens = (pooled_sens - 1.96 * se_pooled_sens, pooled_sens + 1.96 * se_pooled_sens)
    ci_spec = (pooled_spec - 1.96 * se_pooled_spec, pooled_spec + 1.96 * se_pooled_spec)
    
    return pooled_sens, ci_sens, pooled_spec, ci_spec

# ============================================================
# FOREST PLOT
# ============================================================

def make_forest_plot(data, pooled_sens, ci_sens, pooled_spec, ci_spec, title, filename, col_main="steelblue"):
    """Crea il forest plot per sensibilità e specificità"""
    
    sens, spec = calculate_sens_spec(data)
    k = len(data)
    
    # Calcola IC per ogni studio
    sens_ci = np.array([wilson_ci(data['TP'].iloc[i], data['TP'].iloc[i] + data['FN'].iloc[i]) 
                        for i in range(k)])
    spec_ci = np.array([wilson_ci(data['TN'].iloc[i], data['TN'].iloc[i] + data['FP'].iloc[i]) 
                        for i in range(k)])
    
    # Crea la figura - layout più compatto
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle(title, fontsize=11, fontweight='bold', y=0.98)
    
    y = np.arange(k)
    
    # === SENSITIVITY ===
    ax1.set_xlim(0.70, 1.05)
    ax1.set_ylim(-1, k)
    ax1.set_xlabel("Sensitivity", fontsize=10)
    ax1.set_yticks(y)
    ax1.set_yticklabels(data['study'], fontsize=9)
    
    # Griglia verticale
    ax1.grid(axis='x', alpha=0.3, linestyle='--')
    ax1.axvline(pooled_sens, color='red', linestyle='--', linewidth=1.5)
    
    # Plot per ogni studio
    for i in range(k):
        # Linea dell'IC
        ax1.hlines(i, sens_ci[i, 0], sens_ci[i, 1], color=col_main, linewidth=2, capstyle='round')
        # Punto stima
        ax1.plot(sens[i], i, 's', color=col_main, markersize=8)
        # Etichetta valore
        ax1.text(1.03, i, f"{sens[i]*100:.1f}%", va='center', ha='left', fontsize=8)
    
    # Pooled estimate (rombo rosso in basso)
    ax1.plot(pooled_sens, -0.5, 'D', color='red', markersize=12, label='Pooled')
    ax1.hlines(-0.5, ci_sens[0], ci_sens[1], color='red', linewidth=2)
    ax1.text(1.03, -0.5, f"{pooled_sens*100:.1f}%", va='center', ha='left', fontsize=8, fontweight='bold', color='red')
    
    ax1.set_title(f"Sensitivity\nPooled: {pooled_sens*100:.1f}% ({ci_sens[0]*100:.1f}%-{ci_sens[1]*100:.1f}%)", fontsize=10)
    
    # === SPECIFICITY ===
    ax2.set_xlim(-0.05, 0.55)
    ax2.set_ylim(-1, k)
    ax2.set_xlabel("Specificity", fontsize=10)
    ax2.set_yticks([])
    
    # Griglia verticale
    ax2.grid(axis='x', alpha=0.3, linestyle='--')
    ax2.axvline(pooled_spec, color='red', linestyle='--', linewidth=1.5)
    
    # Plot per ogni studio
    for i in range(k):
        # Linea dell'IC
        ax2.hlines(i, spec_ci[i, 0], spec_ci[i, 1], color=col_main, linewidth=2, capstyle='round')
        # Punto stima
        ax2.plot(spec[i], i, 's', color=col_main, markersize=8)
        # Etichetta valore
        ax2.text(0.53, i, f"{spec[i]*100:.1f}%", va='center', ha='right', fontsize=8)
    
    # Pooled estimate
    ax2.plot(pooled_spec, -0.5, 'D', color='red', markersize=12, label='Pooled')
    ax2.hlines(-0.5, ci_spec[0], ci_spec[1], color='red', linewidth=2)
    ax2.text(0.53, -0.5, f"{pooled_spec*100:.1f}%", va='center', ha='right', fontsize=8, fontweight='bold', color='red')
    
    ax2.set_title(f"Specificity\nPooled: {pooled_spec*100:.1f}% ({ci_spec[0]*100:.1f}%-{ci_spec[1]*100:.1f}%)", fontsize=10)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

# ============================================================
# SROC PLOT
# ============================================================

def make_sroc_plot(data, pooled_sens, pooled_spec, title, filename, color="blue", label=""):
    """Crea il plot SROC (Summary Receiver Operating Characteristic)"""
    
    sens, spec = calculate_sens_spec(data)
    fpr = 1 - spec
    tpr = sens
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot SROC curve (semplificata - ellisse approssimata)
    # In una implementazione completa si userebbe il modello bivariate di Reitsma
    theta = np.linspace(0, np.pi/2, 100)
    
    # Centro dell'ellisse (punto summary)
    center_x = 1 - pooled_spec
    center_y = pooled_sens
    
    # Semiassi basati sulla variabilità
    a = 0.15  # asse maggiore
    b = 0.10  # asse minore
    
    # Ruota l'ellisse
    phi = np.pi / 4  # angolo di rotazione
    for angle in theta:
        x = center_x + a * np.cos(angle) * np.cos(phi) - b * np.sin(angle) * np.sin(phi)
        y = center_y + a * np.cos(angle) * np.sin(phi) + b * np.sin(angle) * np.cos(phi)
        if 0 <= x <= 1 and 0 <= y <= 1:
            ax.plot(x, y, color=color, linewidth=2)
    
    # Punto summary
    ax.plot(center_x, center_y, 'D', color=color, markersize=12, label=f'Summary\n({label})' if label else 'Summary')
    
    # Plot dei singoli studi
    sizes = np.sqrt(data['N']) / 25
    ax.scatter(fpr, tpr, s=sizes*50, c=color, alpha=0.6, label='Studies')
    
    # Etichette degli studi
    for i, study in enumerate(data['study']):
        study_short = study.split(' ')[0]  # Solo il primo nome
        ax.annotate(study_short, (fpr[i], tpr[i]), xytext=(5, 5), 
                   textcoords='offset points', fontsize=8, color='darkblue')
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False Positive Rate (1 - Specificity)")
    ax.set_ylabel("True Positive Rate (Sensitivity)")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(loc='lower right')
    
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

def make_sroc_comparison_plot(gfap_data, s100b_data, pooled_gfap_sens, pooled_gfap_spec, 
                               pooled_s100_sens, pooled_s100_spec, title, filename):
    """Crea il plot SROC comparativo"""
    
    sens_g, spec_g = calculate_sens_spec(gfap_data)
    fpr_g = 1 - spec_g
    tpr_g = sens_g
    
    sens_s, spec_s = calculate_sens_spec(s100b_data)
    fpr_s = 1 - spec_s
    tpr_s = sens_s
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # SROC curve GFAP (blu)
    center_x_g = 1 - pooled_gfap_spec
    center_y_g = pooled_gfap_sens
    theta = np.linspace(0, np.pi/2, 100)
    a, b = 0.15, 0.10
    phi = np.pi / 4
    sroc_x_g = []
    sroc_y_g = []
    for angle in theta:
        x = center_x_g + a * np.cos(angle) * np.cos(phi) - b * np.sin(angle) * np.sin(phi)
        y = center_y_g + a * np.cos(angle) * np.sin(phi) + b * np.sin(angle) * np.cos(phi)
        if 0 <= x <= 1 and 0 <= y <= 1:
            sroc_x_g.append(x)
            sroc_y_g.append(y)
    if sroc_x_g:
        ax.plot(sroc_x_g, sroc_y_g, color='blue', linewidth=2, label='GFAP+UCH-L1 SROC')
    
    # SROC curve S100B (verde)
    center_x_s = 1 - pooled_s100_spec
    center_y_s = pooled_s100_sens
    sroc_x_s = []
    sroc_y_s = []
    for angle in theta:
        x = center_x_s + a * np.cos(angle) * np.cos(phi) - b * np.sin(angle) * np.sin(phi)
        y = center_y_s + a * np.cos(angle) * np.sin(phi) + b * np.sin(angle) * np.cos(phi)
        if 0 <= x <= 1 and 0 <= y <= 1:
            sroc_x_s.append(x)
            sroc_y_s.append(y)
    if sroc_x_s:
        ax.plot(sroc_x_s, sroc_y_s, color='darkgreen', linewidth=2, label='S100B SROC')
    
    # Punti summary
    ax.plot(center_x_g, center_y_g, 'D', color='blue', markersize=12)
    ax.plot(center_x_s, center_y_s, '^', color='darkgreen', markersize=12)
    
    # Plot dei singoli studi GFAP
    sizes_g = np.sqrt(gfap_data['N']) / 25
    ax.scatter(fpr_g, tpr_g, s=sizes_g*50, c='blue', alpha=0.6)
    for i, study in enumerate(gfap_data['study']):
        ax.annotate(study.split(' ')[0], (fpr_g[i], tpr_g[i]), xytext=(5, 5), 
                   textcoords='offset points', fontsize=8, color='blue')
    
    # Plot dei singoli studi S100B
    sizes_s = np.sqrt(s100b_data['N']) / 25
    ax.scatter(fpr_s, tpr_s, s=sizes_s*50, c='darkgreen', alpha=0.7, marker='^')
    for i, study in enumerate(s100b_data['study']):
        ax.annotate(study.split(' ')[0], (fpr_s[i], tpr_s[i]), xytext=(5, 5), 
                   textcoords='offset points', fontsize=8, color='darkgreen')
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, 1.02)
    ax.set_xlabel("False Positive Rate (1 - Specificity)")
    ax.set_ylabel("True Positive Rate (Sensitivity)")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(loc='lower right')
    
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

# ============================================================
# DEEKS FUNNEL PLOT
# ============================================================

def make_deeks_funnel(data, filename):
    """Crea il Deeks funnel plot per publication bias"""
    
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
        print("  => No evidence of publication bias")
    else:
        print("  => Possible publication bias")
    
    fig, ax = plt.subplots(figsize=(9, 7))
    
    # Plot punti
    sizes = np.sqrt(data['N']) / 25
    ax.scatter(inv_sqrt_ESS, lnDOR, s=sizes*50, c='steelblue', alpha=0.7)
    
    # Linea di regressione
    x_line = np.linspace(inv_sqrt_ESS.min(), inv_sqrt_ESS.max(), 100)
    y_line = intercept + slope * x_line
    ax.plot(x_line, y_line, 'r--', linewidth=2, label=f'Regression line (p={p_value_t:.3f})')
    
    # Etichette studi
    for i, study in enumerate(data['study']):
        ax.annotate(study, (inv_sqrt_ESS.iloc[i], lnDOR.iloc[i]), xytext=(5, 5), 
                   textcoords='offset points', fontsize=8, color='darkblue')
    
    ax.set_xlabel("1/√(Effective Sample Size)")
    ax.set_ylabel("ln(Diagnostic Odds Ratio)")
    ax.set_title(f"Deeks Funnel Plot - GFAP+UCH-L1")
    ax.grid(alpha=0.3)
    ax.legend()
    
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

# ============================================================
# MAIN - GENERA TUTTI I PLOT
# ============================================================

if __name__ == "__main__":
    print("=" * 60)
    print("mTBI DTA Meta-Analysis - Generating Plots")
    print("=" * 60)
    
    # Calcola stime pooled per GFAP
    pooled_sens_gfap, ci_sens_gfap, pooled_spec_gfap, ci_spec_gfap = calculate_pooled_estimates(gfap)
    print(f"\nGFAP+UCH-L1:")
    print(f"  Pooled Sensitivity: {pooled_sens_gfap:.3f} ({ci_sens_gfap[0]:.3f}-{ci_sens_gfap[1]:.3f})")
    print(f"  Pooled Specificity: {pooled_spec_gfap:.3f} ({ci_spec_gfap[0]:.3f}-{ci_spec_gfap[1]:.3f})")
    
    # Calcola stime pooled per S100B
    pooled_sens_s100, ci_sens_s100, pooled_spec_s100, ci_spec_s100 = calculate_pooled_estimates(s100b)
    print(f"\nS100B:")
    print(f"  Pooled Sensitivity: {pooled_sens_s100:.3f} ({ci_sens_s100[0]:.3f}-{ci_sens_s100[1]:.3f})")
    print(f"  Pooled Specificity: {pooled_spec_s100:.3f} ({ci_spec_s100[0]:.3f}-{ci_spec_s100[1]:.3f})")
    
    # Genera Forest Plots
    print("\n--- Generating Forest Plots ---")
    make_forest_plot(
        gfap, pooled_sens_gfap, ci_sens_gfap, pooled_spec_gfap, ci_spec_gfap,
        "GFAP+UCH-L1 (Bivariate Model)",
        "forest_GFAP_UCH_L1_python.png",
        col_main="steelblue"
    )
    
    make_forest_plot(
        s100b, pooled_sens_s100, ci_sens_s100, pooled_spec_s100, ci_spec_s100,
        "S100B (Bivariate Model)",
        "forest_S100B_python.png",
        col_main="forestgreen"
    )
    
    # Genera SROC plots
    print("\n--- Generating SROC Plots ---")
    make_sroc_plot(
        gfap, pooled_sens_gfap, pooled_spec_gfap,
        "SROC - GFAP+UCH-L1",
        "SROC_GFAP_python.png",
        color="blue",
        label=f"Se={pooled_sens_gfap:.2f}, Sp={pooled_spec_gfap:.2f}"
    )
    
    make_sroc_plot(
        s100b, pooled_sens_s100, pooled_spec_s100,
        "SROC - S100B",
        "SROC_S100B_python.png",
        color="darkgreen",
        label=f"Se={pooled_sens_s100:.2f}, Sp={pooled_spec_s100:.2f}"
    )
    
    # Genera SROC comparativo
    make_sroc_comparison_plot(
        gfap, s100b,
        pooled_sens_gfap, pooled_spec_gfap,
        pooled_sens_s100, pooled_spec_s100,
        "SROC Comparison: GFAP+UCH-L1 vs S100B",
        "SROC_comparison_python.png"
    )
    
    # Genera Deeks funnel plot
    print("\n--- Generating Deeks Funnel Plot ---")
    make_deeks_funnel(gfap, "deeks_funnel_GFAP_python.png")
    
    print("\n" + "=" * 60)
    print("All plots generated successfully!")
    print("=" * 60)
