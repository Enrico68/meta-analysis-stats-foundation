# ============================================================
# mTBI Biomarkers DTA Meta-Analysis — v5 (FIXED)
# Fix definitivo: uso esplicito mada:: per evitare conflitti con metafor
# ============================================================

if (!require("mada")) install.packages("mada", repos="https://cran.r-project.org")
if (!require("metafor")) install.packages("metafor", repos="https://cran.r-project.org")

library(mada)      # carichiamo solo mada (metafor lo useremo solo per Deeks)
# ============================================================
# 1. DATA
# ============================================================
gfap <- data.frame(
  study = c("Bazarian 2021", "Lagares 2024", "Faisal 2023",
            "Legramante 2024", "Ladang 2025", "Milevoj K. 2025",
            "Curran 2025", "Campagna 2025"),
  TP = c(115, 176, 12, 7, 112, 107, 12, 8),
  FP = c(1061, 946, 75, 89, 148, 529, 56, 59),
  FN = c(5, 3, 0, 0, 1, 5, 0, 0),
  TN = c(720, 313, 43, 34, 101, 181, 21, 30)
)
gfap$N <- gfap$TP + gfap$FP + gfap$FN + gfap$TN

s100b <- data.frame(
  study = c("Oris 2024", "Oris 2021", "Hopman 2023"),
  TP = c(7, 63, 69),
  FP = c(148, 891, 355),
  FN = c(0, 0, 5),
  TN = c(51, 218, 66)
)
s100b$N <- s100b$TP + s100b$FP + s100b$FN + s100b$TN

gfap_sens <- data.frame(
  study = c("Lapic 2025", "Jalali 2025", "Spaziani 2025", "Welch 2025"),
  TP = c(27, 5, 72, 116),
  FP = c(104, 83, 277, 1066),
  FN = c(2, 1, 3, 4),
  TN = c(42, 34, 89, 713)
)
gfap_sens$N <- gfap_sens$TP + gfap_sens$FP + gfap_sens$FN + gfap_sens$TN

# ============================================================
# 2. REITSMA MODELS
# ============================================================
cat("\n==============================================================\n")
cat(" GFAP+UCH-L1 — BIVARIATE REITSMA (8 studies, N =", sum(gfap$N), ")\n")
cat("==============================================================\n\n")

fit_gfap <- mada::reitsma(gfap, formula = cbind(tsens, tfpr) ~ 1,
                          correction.control = "all", correction = 0.5)

print(summary(fit_gfap))

cat("\n--- POOLED ESTIMATES ---\n")
cat(sprintf("Sensitivity: %.3f (95%% CI: %.3f - %.3f)\n",
            fit_gfap$Sens, fit_gfap$Sens.ci[1], fit_gfap$Sens.ci[2]))
cat(sprintf("Specificity: %.3f (95%% CI: %.3f - %.3f)\n",
            fit_gfap$Spec, fit_gfap$Spec.ci[1], fit_gfap$Spec.ci[2]))
cat(sprintf("AUC: %.3f | pAUC: %.3f\n", summary(fit_gfap)$AUC, summary(fit_gfap)$pAUC))

npv <- sum(gfap$TN) / (sum(gfap$TN) + sum(gfap$FN))
ct_red <- (sum(gfap$FN) + sum(gfap$TN)) / sum(gfap$N) * 100
cat(sprintf("NPV: %.3f | CT avoided: %.1f%% | Missed: %d/%d (%.1f%%)\n",
            npv, ct_red, sum(gfap$FN), sum(gfap$TP+gfap$FN),
            sum(gfap$FN)/sum(gfap$TP+gfap$FN)*100))

# ============================================================
# 3. FOREST PLOTS (usa esplicitamente mada::)
# ============================================================
cat("\n\n--- Generating Plots ---\n")

png("forest_GFAP_UCH_L1_reitsma.png", width=1400, height=800, res=150)
mada::forest(fit_gfap,
             main = "GFAP+UCH-L1: Paired Forest Plot (Reitsma)",
             snames = gfap$study)
dev.off()
cat("Saved: forest_GFAP_UCH_L1_reitsma.png\n")

# ============================================================
# 4. SROC CURVES
# ============================================================
png("SROC_GFAP_reitsma.png", width=900, height=900, res=150)
mada::plot(fit_gfap, sroclwd=2, predict=TRUE, conf.region=TRUE,
           main="SROC - GFAP+UCH-L1 (Bivariate Reitsma)")
fpr_g <- 1 - gfap$TN/(gfap$TN+gfap$FP)
tpr_g <- gfap$TP/(gfap$TP+gfap$FN + 0.5*(gfap$FN==0))
points(fpr_g, tpr_g, pch=19, cex=sqrt(gfap$N)/25, col=rgb(0.1,0.3,0.8,0.6))
text(fpr_g, tpr_g, labels=gsub(" \\d+$","",gfap$study), pos=4, cex=0.55, col="darkblue")
dev.off()
cat("Saved: SROC_GFAP_reitsma.png\n")

# ============================================================
# 5. DEEKS FUNNEL (ri-carichiamo metafor solo qui)
# ============================================================
library(metafor)
gfap$DOR <- ((gfap$TP+0.5)*(gfap$TN+0.5)) / ((gfap$FP+0.5)*(gfap$FN+0.5))
gfap$lnDOR <- log(gfap$DOR)
gfap$var_lnDOR <- 1/(gfap$TP+0.5)+1/(gfap$FP+0.5)+1/(gfap$FN+0.5)+1/(gfap$TN+0.5)
gfap$ESS <- 4*gfap$N / (((gfap$TP+gfap$FN)/gfap$N)+((gfap$FP+gfap$TN)/gfap$N))^2
gfap$inv_sqrt_ESS <- 1/sqrt(gfap$ESS)

deeks <- lm(lnDOR ~ inv_sqrt_ESS, data=gfap, weights=1/var_lnDOR)
deeks_p <- summary(deeks)$coefficients["inv_sqrt_ESS","Pr(>|t|)"]

cat(sprintf("\nDeeks test p = %.3f", deeks_p))
if (deeks_p > 0.10) cat(" => No publication bias\n") else cat(" => Possible bias\n")

png("deeks_funnel_GFAP.png", width=900, height=700, res=150)
plot(gfap$inv_sqrt_ESS, gfap$lnDOR, pch=19, cex=sqrt(gfap$N)/25, col="steelblue",
     xlab="1/sqrt(ESS)", ylab="ln(DOR)",
     main=sprintf("Deeks Funnel Plot - GFAP+UCH-L1 (p=%.3f)", deeks_p))
abline(deeks, col="red", lwd=2, lty=2)
text(gfap$inv_sqrt_ESS, gfap$lnDOR, labels=gfap$study, pos=4, cex=0.6)
grid()
dev.off()
cat("Saved: deeks_funnel_GFAP.png\n")

# ============================================================
# 6. FINAL SUMMARY FOR MANUSCRIPT
# ============================================================
cat("\n\n==============================================================\n")
cat(" FINAL RESULTS FOR MANUSCRIPT (JAMA Network Open)\n")
cat("==============================================================\n\n")
cat("GFAP+UCH-L1 (8 studies, N=4969):\n")
cat(sprintf(" Se %.1f%% (%.1f-%.1f) | Sp %.1f%% (%.1f-%.1f)\n",
            fit_gfap$Sens*100, fit_gfap$Sens.ci[1]*100, fit_gfap$Sens.ci[2]*100,
            fit_gfap$Spec*100, fit_gfap$Spec.ci[1]*100, fit_gfap$Spec.ci[2]*100))
cat(sprintf(" AUC=%.3f | NPV=%.1f%% | CT reduction=%.1f%%\n",
            summary(fit_gfap)$AUC, npv*100, ct_red))
cat("\n=== ANALYSIS COMPLETE ===\n")