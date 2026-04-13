# =============================================================================
# DATA INPUT — DATASET v15 DEFINITIVO (tutti verificati da full text)
# Data: 11 aprile 2026
# =============================================================================

# --- GFAP+UCH-L1: 7 studi (pool primario inclusivo, tutte le piattaforme) ---
GFAP_data <- data.frame(
  study = c("Chayoua 2024", "Bazarian 2021", "Papa 2022",
            "Lagares 2024 (BRAINI)", "Puravet 2026",
            "Milevoj 2025", "Legramante 2024"),
  TP = c( 57, 115,  23, 176,  50, 107,   7),
  FP = c(157,1061, 245, 945, 700, 529,  89),
  FN = c(  2,   5,   0,   3,   0,   5,   0),
  TN = c( 37, 720,  81, 314, 260, 181,  34)
)
# Note v15:
# - Oris 2024 IJMS RIMOSSO (sottoinsieme Puravet 2026)
# - Jones 2020 RIMOSSO (era S100B, non GFAP)
# - Puravet 2025 AEM RIMOSSO (meta-analisi, non studio primario)
# - Chayoua: TP=57 FP=157 FN=2 TN=37 (corretto da full text, v14 aveva dati errati)
# - Bazarian: N=1901, TP=115 (corretto da full text, v14 aveva N=391)
# - Papa: N=349, cutoff 67/189 ELISA (corretto, v14 aveva N=858)
# - Lagares: N=1438, TP=176 (corretto da full text, v14 aveva N=1210)
# - Milevoj 2025: confermato da full text
# - Legramante 2024: confermato da full text

# --- S100B: 4 studi (pool primario) ---
S100B_data <- data.frame(
  study = c("Puravet 2026", "Seidenfaden 2021 (PreTBI I)",
            "Rogan 2023", "Hopman 2023"),
  TP = c( 48,  32,  15,  69),
  FP = c(720, 366,  81, 355),
  FN = c(  2,   0,   1,   5),
  TN = c(240, 168,  36,  66)
)
# Note v15:
# - Oris 2024 CCLM RIMOSSO (sottoinsieme Puravet 2026)
# - Oris 2021 RIMOSSO (overlap Clermont-Ferrand con Puravet 2026)
# - Jones 2020 RIMOSSO (Se=84.6% outlier, cutoff race-specific) → sensibilità
# - Seidenfaden: dati IN-HOSPITAL (Table 5), non preospedalieri
# - Hopman 2023: dati ricostruiti da abstract (full text non disponibile)
