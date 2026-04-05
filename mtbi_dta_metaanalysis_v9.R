################################################################################
# mTBI Biomarker DTA Meta-Analysis — Version 9
# Blood Biomarkers in the Emergency Department Management of mTBI
# Aggiornamento: Pool S100B definitivo a 6 studi (N=3,582)
# Pool GFAP+UCH-L1: 8 studi (N=4,969) — invariato
# Metodo: Reitsma (mada), forest plots, SROC, clinical summary
# Data: Aprile 2026
################################################################################

# ---- 0. Installazione e caricamento pacchetti --------------------------------
if (!requireNamespace("mada", quietly = TRUE)) install.packages("mada")
if (!requireNamespace("meta", quietly = TRUE)) install.packages("meta")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

library(mada)
library(meta)
library(ggplot2)
library(gridExtra)

################################################################################
# 1. DATASET
################################################################################

# ---- 1a. Pool GFAP+UCH-L1 (8 studi, N=4,969) — INVARIATO -------------------
gfap_data <- data.frame(
  study = c(
    "Oris 2024 (CCLM)",
    "Bazarian 2021 (AEM)",
    "Papa 2022 (JAMA NO)",
    "Chayoua 2024 (JNT)",
    "Legramante 2024 (IJEM)",
    "Jones 2020 (Brain Inj)",
    "Lagares 2024 (EBioMed)",
    "Puravet 2025 (AEM)"
  ),
  TP = c(18, 24, 52, 38, 7, 33, 68, 99),
  FP = c(61, 199, 392, 130, 89, 425, 477, 285),
  FN = c(0, 1, 1, 1, 0, 0, 2, 1),
  TN = c(127, 167, 413, 149, 34, 221, 663, 373),
  stringsAsFactors = FALSE
)

gfap_data$N    <- gfap_data$TP + gfap_data$FP + gfap_data$FN + gfap_data$TN
gfap_data$Se   <- gfap_data$TP / (gfap_data$TP + gfap_data$FN)
gfap_data$Sp   <- gfap_data$TN / (gfap_data$TN + gfap_data$FP)
gfap_data$CT_pos <- gfap_data$TP + gfap_data$FN
gfap_data$prev   <- gfap_data$CT_pos / gfap_data$N

cat("=== GFAP+UCH-L1 POOL (8 studi) ===\n")
cat("N totale:", sum(gfap_data$N), "\n")
cat("CT+ totali:", sum(gfap_data$CT_pos), "\n\n")

# ---- 1b. Pool S100B (6 studi, N=3,582) — AGGIORNATO v9 ----------------------
# Studio                 TP   FP  FN   TN    N
# Oris 2024             7   148   0   51  206
# Oris 2021            63   891   0  218 1172
# Hopman 2023          69   355   5   66  495
# Seidenfaden 2021     32   366   0  168  566
# Rogan 2023           15    81   1   36  133
# Puravet 2025/CCLM    48   720   2  240 1010
s100b_data <- data.frame(
  study = c(
    "Oris 2024 (CCLM)",
    "Oris 2021 (CCLM)",
    "Hopman 2023 (Brain Inj)",
    "Seidenfaden 2021 (PreTBI I)",
    "Rogan 2023 (BRAIN, EMJ)",
    "Puravet 2026 (CCLM)"
  ),
  TP = c(7,  63, 69, 32, 15, 48),
  FP = c(148, 891, 355, 366, 81, 720),
  FN = c(0,   0,  5,  0,  1,  2),
  TN = c(51, 218, 66, 168, 36, 240),
  stringsAsFactors = FALSE
)

s100b_data$N    <- s100b_data$TP + s100b_data$FP + s100b_data$FN + s100b_data$TN
s100b_data$Se   <- s100b_data$TP / (s100b_data$TP + s100b_data$FN)
s100b_data$Sp   <- s100b_data$TN / (s100b_data$TN + s100b_data$FP)
s100b_data$CT_pos <- s100b_data$TP + s100b_data$FN
s100b_data$prev   <- s100b_data$CT_pos / s100b_data$N

cat("=== S100B POOL (6 studi) ===\n")
cat("N totale:", sum(s100b_data$N), "\n")
cat("CT+ totali:", sum(s100b_data$CT_pos), "\n\n")

################################################################################
# 2. MODELLO REITSMA (BIVARIATE)
################################################################################

cat("=====================================\n")
cat("MODELLO REITSMA — GFAP+UCH-L1\n")
cat("=====================================\n")

fit_gfap <- reitsma(
  TP = gfap_data$TP, FP = gfap_data$FP,
  FN = gfap_data$FN, TN = gfap_data$TN,
  data = gfap_data
)
print(summary(fit_gfap))

# Stime pooled GFAP
se_gfap <- plogis(fit_gfap$coefficients[1])
sp_gfap <- plogis(fit_gfap$coefficients[2])
cat("\nSensibilità pooled GFAP:", round(se_gfap * 100, 1), "%\n")
cat("Specificità pooled GFAP:", round(sp_gfap * 100, 1), "%\n")

cat("\n=====================================\n")
cat("MODELLO REITSMA — S100B\n")
cat("=====================================\n")

fit_s100b <- reitsma(
  TP = s100b_data$TP, FP = s100b_data$FP,
  FN = s100b_data$FN, TN = s100b_data$TN,
  data = s100b_data
)
print(summary(fit_s100b))

# Stime pooled S100B
se_s100b <- plogis(fit_s100b$coefficients[1])
sp_s100b <- plogis(fit_s100b$coefficients[2])
cat("\nSensibilità pooled S100B:", round(se_s100b * 100, 1), "%\n")
cat("Specificità pooled S100B:", round(sp_s100b * 100, 1), "%\n")

################################################################################
# 3. NPV E CLINICAL METRICS (PREVALENZA MEDIA PONDERATA)
################################################################################

# Prevalenza CT+ ponderata per N
prev_gfap_weighted  <- sum(gfap_data$CT_pos)  / sum(gfap_data$N)
prev_s100b_weighted <- sum(s100b_data$CT_pos) / sum(s100b_data$N)

cat("\n=== PREVALENZE CT+ PONDERATE ===\n")
cat("GFAP+UCH-L1:", round(prev_gfap_weighted * 100, 1), "%\n")
cat("S100B:       ", round(prev_s100b_weighted * 100, 1), "%\n\n")

# NPV con formula di Bayes
npv <- function(se, sp, prev) {
  (sp * (1 - prev)) / (sp * (1 - prev) + (1 - se) * prev)
}

npv_gfap  <- npv(se_gfap,  sp_gfap,  prev_gfap_weighted)
npv_s100b <- npv(se_s100b, sp_s100b, prev_s100b_weighted)

cat("NPV GFAP+UCH-L1:", round(npv_gfap * 100, 2), "%\n")
cat("NPV S100B:       ", round(npv_s100b * 100, 2), "%\n\n")

# Lesioni mancate per 1000 pazienti
miss_gfap  <- (1 - se_gfap)  * prev_gfap_weighted  * 1000
miss_s100b <- (1 - se_s100b) * prev_s100b_weighted * 1000

cat("Lesioni mancate/1000 pz — GFAP+UCH-L1:", round(miss_gfap, 1), "\n")
cat("Lesioni mancate/1000 pz — S100B:       ", round(miss_s100b, 1), "\n\n")

# CT evitate (pazienti con test negativo)
ct_avoid_gfap  <- sp_gfap  * (1 - prev_gfap_weighted)  * 100
ct_avoid_s100b <- sp_s100b * (1 - prev_s100b_weighted) * 100

cat("CT potenzialmente evitate — GFAP+UCH-L1:", round(ct_avoid_gfap, 1), "%\n")
cat("CT potenzialmente evitate — S100B:       ", round(ct_avoid_s100b, 1), "%\n\n")

################################################################################
# 4. FOREST PLOT — GFAP+UCH-L1
################################################################################

png("Fig2_forest_gfap_v9.png", width = 2800, height = 1800, res = 200)

par(mfrow = c(1, 2), mar = c(5, 12, 4, 2), oma = c(0, 0, 3, 0))

# Sensibilità
forest(
  madad(TP = gfap_data$TP, FP = gfap_data$FP,
        FN = gfap_data$FN, TN = gfap_data$TN),
  type = "sens",
  main = "Sensibilità",
  xlab = "Sensibilità (95% CI)",
  snames = gfap_data$study,
  xlim = c(0.5, 1.05),
  col.diamond = "#2563EB",
  col.square  = "#93C5FD"
)
abline(v = se_gfap, col = "#2563EB", lty = 2, lwd = 2)

# Specificità
forest(
  madad(TP = gfap_data$TP, FP = gfap_data$FP,
        FN = gfap_data$FN, TN = gfap_data$TN),
  type = "spec",
  main = "Specificità",
  xlab = "Specificità (95% CI)",
  snames = gfap_data$study,
  xlim = c(0, 0.7),
  col.diamond = "#2563EB",
  col.square  = "#93C5FD"
)
abline(v = sp_gfap, col = "#2563EB", lty = 2, lwd = 2)

mtext("Forest Plot — GFAP+UCH-L1 (8 studi, N=4,969)",
      outer = TRUE, cex = 1.4, font = 2, col = "#1e3a5f")
dev.off()
cat("Salvato: Fig2_forest_gfap_v9.png\n")

################################################################################
# 5. FOREST PLOT — S100B
################################################################################

png("Fig3_forest_s100b_v9.png", width = 2800, height = 1900, res = 200)

par(mfrow = c(1, 2), mar = c(5, 15, 4, 2), oma = c(0, 0, 3, 0))

# Sensibilità
forest(
  madad(TP = s100b_data$TP, FP = s100b_data$FP,
        FN = s100b_data$FN, TN = s100b_data$TN),
  type = "sens",
  main = "Sensibilità",
  xlab = "Sensibilità (95% CI)",
  snames = s100b_data$study,
  xlim = c(0.5, 1.1),
  col.diamond = "#16A34A",
  col.square  = "#86EFAC"
)
abline(v = se_s100b, col = "#16A34A", lty = 2, lwd = 2)

# Specificità
forest(
  madad(TP = s100b_data$TP, FP = s100b_data$FP,
        FN = s100b_data$FN, TN = s100b_data$TN),
  type = "spec",
  main = "Specificità",
  xlab = "Specificità (95% CI)",
  snames = s100b_data$study,
  xlim = c(0, 0.5),
  col.diamond = "#16A34A",
  col.square  = "#86EFAC"
)
abline(v = sp_s100b, col = "#16A34A", lty = 2, lwd = 2)

mtext("Forest Plot — S100B (6 studi, N=3,582)",
      outer = TRUE, cex = 1.4, font = 2, col = "#14532d")
dev.off()
cat("Salvato: Fig3_forest_s100b_v9.png\n")

################################################################################
# 6. SROC COMPARATIVO
################################################################################

png("Fig4_sroc_v9.png", width = 2400, height = 2200, res = 200)

par(mar = c(5, 5, 4, 2))

# SROC GFAP
plot(fit_gfap,
     xlim = c(0, 1), ylim = c(0, 1),
     main = "SROC — GFAP+UCH-L1 vs S100B",
     xlab = "1 − Specificità (Tasso Falsi Positivi)",
     ylab = "Sensibilità (Tasso Veri Positivi)",
     col  = "#2563EB",
     lwd  = 2.5,
     cex.main = 1.3,
     cex.lab  = 1.1
)

# SROC S100B (aggiunto)
lines(sroc(fit_s100b), col = "#16A34A", lwd = 2.5)

# Punti sommario
points(1 - sp_gfap,  se_gfap,  pch = 23, bg = "#2563EB", cex = 2.5, col = "white")
points(1 - sp_s100b, se_s100b, pch = 22, bg = "#16A34A", cex = 2.5, col = "white")

# Punti individuali
points(1 - gfap_data$Sp,  gfap_data$Se,  pch = 21,
       bg = "#BFDBFE", col = "#2563EB", cex = sqrt(gfap_data$N) / 8)
points(1 - s100b_data$Sp, s100b_data$Se, pch = 21,
       bg = "#BBF7D0", col = "#16A34A", cex = sqrt(s100b_data$N) / 8)

# Regione CI ellittica
lines(fit_gfap,  "conf",  col = "#2563EB", lty = 2)
lines(fit_s100b, "conf",  col = "#16A34A", lty = 2)

# Legenda
legend("bottomright",
       legend = c(
         sprintf("GFAP+UCH-L1: Se=%.1f%%, Sp=%.1f%%", se_gfap*100, sp_gfap*100),
         sprintf("S100B:        Se=%.1f%%, Sp=%.1f%%", se_s100b*100, sp_s100b*100)
       ),
       col = c("#2563EB", "#16A34A"),
       lwd = c(2.5, 2.5),
       pch = c(23, 22),
       pt.bg = c("#2563EB", "#16A34A"),
       pt.cex = 1.8,
       bty = "n",
       cex = 1.0
)

abline(a = 0, b = 1, lty = 3, col = "grey60")
dev.off()
cat("Salvato: Fig4_sroc_v9.png\n")

################################################################################
# 7. CLINICAL SUMMARY FIGURE (ggplot2)
################################################################################

# Tabella riepilogativa
summary_df <- data.frame(
  Biomarker    = c("GFAP+UCH-L1", "S100B"),
  N_studi      = c(nrow(gfap_data), nrow(s100b_data)),
  N_pazienti   = c(sum(gfap_data$N), sum(s100b_data$N)),
  Sensibilita  = c(round(se_gfap  * 100, 1), round(se_s100b  * 100, 1)),
  Specificita  = c(round(sp_gfap  * 100, 1), round(sp_s100b  * 100, 1)),
  NPV          = c(round(npv_gfap  * 100, 2), round(npv_s100b  * 100, 2)),
  CT_evitate   = c(round(ct_avoid_gfap, 1),  round(ct_avoid_s100b, 1)),
  Miss_per1000 = c(round(miss_gfap, 1),       round(miss_s100b, 1))
)

cat("\n=== TABELLA RIEPILOGATIVA FINALE ===\n")
print(summary_df)

# Plot comparativo a barre
metrics_long <- data.frame(
  Biomarker = rep(c("GFAP+UCH-L1", "S100B"), each = 4),
  Metrica   = rep(c("Sensibilità (%)", "Specificità (%)",
                    "NPV (%)", "CT evitate (%)"), 2),
  Valore    = c(
    se_gfap*100,  sp_gfap*100,  npv_gfap*100,  ct_avoid_gfap,
    se_s100b*100, sp_s100b*100, npv_s100b*100, ct_avoid_s100b
  )
)

png("Fig5_clinical_summary_v9.png", width = 2800, height = 1800, res = 200)

p <- ggplot(metrics_long, aes(x = Metrica, y = Valore, fill = Biomarker)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, alpha = 0.9) +
  geom_text(aes(label = round(Valore, 1)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("GFAP+UCH-L1" = "#2563EB", "S100B" = "#16A34A")) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(
    title    = "Confronto Clinico: GFAP+UCH-L1 vs S100B nel mTBI",
    subtitle = sprintf(
      "GFAP+UCH-L1: %d studi, N=%d | S100B: %d studi, N=%d",
      nrow(gfap_data), sum(gfap_data$N),
      nrow(s100b_data), sum(s100b_data$N)
    ),
    x     = NULL,
    y     = "Valore (%)",
    fill  = "Biomarker",
    caption = "Modello Reitsma (bivariate random-effects). CT evitate = pazienti con test negativo tra CT−."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    legend.position = "top",
    legend.title  = element_text(face = "bold"),
    axis.text.x   = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    plot.caption  = element_text(size = 9, color = "grey50", hjust = 0)
  )

print(p)
dev.off()
cat("Salvato: Fig5_clinical_summary_v9.png\n")

################################################################################
# 8. SENSITIVITY ANALYSIS — S100B con e senza Puravet 2026
################################################################################

cat("\n=== SENSITIVITY ANALYSIS S100B ===\n")
cat("--- Pool completo (6 studi, incluso Puravet 2026) ---\n")
cat("Se:", round(se_s100b*100,1), "%, Sp:", round(sp_s100b*100,1), "%\n\n")

# Pool senza Puravet (5 studi originali)
s100b_no_pur <- s100b_data[-6, ]
fit_s100b_np <- reitsma(
  TP = s100b_no_pur$TP, FP = s100b_no_pur$FP,
  FN = s100b_no_pur$FN, TN = s100b_no_pur$TN,
  data = s100b_no_pur
)
se_np <- plogis(fit_s100b_np$coefficients[1])
sp_np <- plogis(fit_s100b_np$coefficients[2])
cat("--- Pool senza Puravet (5 studi, N=2,572) ---\n")
cat("Se:", round(se_np*100,1), "%, Sp:", round(sp_np*100,1), "%\n\n")

cat("Differenza Se:", round((se_s100b - se_np)*100, 1), "pp\n")
cat("Differenza Sp:", round((sp_s100b - sp_np)*100, 1), "pp\n")

################################################################################
# 9. OUTPUT FINALE CONSOLEBLE
################################################################################

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║           RISULTATI FINALI — Meta-analisi mTBI v9           ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  GFAP+UCH-L1  │ 8 studi │ N=%d │ Se=%.1f%% │ Sp=%.1f%%  ║\n",
            sum(gfap_data$N), se_gfap*100, sp_gfap*100))
cat(sprintf("║               │         │        │ NPV=%.2f%%            ║\n",
            npv_gfap*100))
cat(sprintf("║               │         │        │ CT evitate: %.1f%%    ║\n",
            ct_avoid_gfap))
cat(sprintf("║               │         │        │ Miss/1000: %.1f       ║\n",
            miss_gfap))
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  S100B        │ 6 studi │ N=%d │ Se=%.1f%% │ Sp=%.1f%%  ║\n",
            sum(s100b_data$N), se_s100b*100, sp_s100b*100))
cat(sprintf("║               │         │        │ NPV=%.2f%%            ║\n",
            npv_s100b*100))
cat(sprintf("║               │         │        │ CT evitate: %.1f%%    ║\n",
            ct_avoid_s100b))
cat(sprintf("║               │         │        │ Miss/1000: %.1f       ║\n",
            miss_s100b))
cat("╚══════════════════════════════════════════════════════════════╝\n")

cat("\nScript completato. File generati:\n")
cat("  - Fig2_forest_gfap_v9.png\n")
cat("  - Fig3_forest_s100b_v9.png\n")
cat("  - Fig4_sroc_v9.png\n")
cat("  - Fig5_clinical_summary_v9.png\n")
cat("\nVersione: v9 | S100B pool: 6 studi, N=3,582\n")
