# 📊 Meta-Analysis Stats Foundation

> **Fondamenti statistici per la Meta-Analisi** — Un corso progressivo in R Notebooks

## 🎯 Obiettivo

Questo repository contiene i notebook fondazionali per apprendere i concetti statistici
necessari per condurre e interpretare una meta-analisi in R.

Progettato per essere usato con **Positron IDE** + **Claude Code**.

---

## 📚 Struttura del Corso

| # | Notebook | Argomento | Stato |
|---|----------|-----------|-------|
| 01 | `01_distribuzioni_e_probabilita.qmd` | Distribuzioni, probabilità, normale | ✅ |
| 02 | `02_stima_e_inferenza.qmd` | Media, varianza, SE, IC | ✅ |
| 03 | `03_effect_size.qmd` | Cohen's d, OR, RR, correlazione | ✅ |
| 04 | `04_eterogeneita.qmd` | Q, I², τ² — cuore della meta-analisi | ✅ |
| 05 | `05_modelli_fixed_random.qmd` | Fixed vs Random effects | ✅ |

---

## 🛠️ Setup

```r
install.packages(c("tidyverse", "meta", "metafor", "ggplot2", "knitr", "effectsize"))
```

## 🚀 Come Usare

1. Apri il progetto in **Positron**
2. Esegui i notebook in ordine (01 → 05)
3. Usa **Claude Code** per esplorare, modificare, sperimentare

---

*Testi di riferimento: Borenstein et al. (2021) · Harrer et al. (2022)*
