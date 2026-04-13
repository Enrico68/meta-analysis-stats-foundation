#!/usr/bin/env Rscript
# =============================================================================
# mTBI Biomarker DTA Meta-Analysis - Version 15 DEFINITIVO (metafor)
# Dataset verificato dai full text PDF — 11 aprile 2026
# Skeleton identico a v14 funzionante, dati aggiornati
# =============================================================================

packages <- c("metafor", "ggplot2", "gridExtra", "grid", "mada")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# =============================================================================
# DATA INPUT — DATASET v15 DEFINITIVO
# =============================================================================

GFAP_data <- data.frame(
  study = c("Chayoua 2024", "Bazarian 2021", "Papa 2022",
            "Lagares 2024 (BRAINI)", "Puravet 2026",
            "Milevoj 2025", "Legramante 2024"),
  TP = c( 57, 115,  23, 176,  50, 107,   7),
  FP = c(157,1061, 245, 945, 700, 529,  89),
  FN = c(  2,   5,   0,   3,   0,   5,   0),
  TN = c( 37, 720,  81, 314, 260, 181,  34)
)

S100B_data <- data.frame(
  study = c("Puravet 2026", "Seidenfaden 2021",
            "Rogan 2023", "Hopman 2023"),
  TP = c( 48,  32,  15,  69),
  FP = c(720, 366,  81, 355),
  FN = c(  2,   0,   1,   5),
  TN = c(240, 168,  36,  66)
)

apply_cc <- function(data, add = 0.5) {
  data$TP <- ifelse(data$TP == 0, data$TP + add, data$TP)
  data$FN <- ifelse(data$FN == 0, data$FN + add, data$FN)
  data$FP <- ifelse(data$FP == 0, data$FP + add, data$FP)
  data$TN <- ifelse(data$TN == 0, data$TN + add, data$TN)
  return(data)
}

GFAP_data  <- apply_cc(GFAP_data)
S100B_data <- apply_cc(S100B_data)

wilson_ci <- function(p, n, z = 1.96) {
  if (n == 0) return(c(lo = 0, hi = 1))
  center <- (p + z^2/(2*n)) / (1 + z^2/n)
  margin <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n)
  return(c(lo = max(0, center - margin), hi = min(1, center + margin)))
}

run_meta_analysis <- function(data, label) {
  data$Se    <- data$TP / (data$TP + data$FN)
  data$Sp    <- data$TN / (data$TN + data$FP)
  data$N     <- data$TP + data$FP + data$FN + data$TN
  data$CT_pos<- data$TP + data$FN
  data$prev  <- data$CT_pos / data$N

  dat_se <- escalc(measure = "PLO", xi = data$TP,
                   ni = data$TP + data$FN, data = data, slab = study)
  fit_se <- rma(yi, vi, data = dat_se, method = "REML")
  pred_se<- predict(fit_se, transf = transf.ilogit)

  dat_sp <- escalc(measure = "PLO", xi = data$TN,
                   ni = data$TN + data$FP, data = data, slab = study)
  fit_sp <- rma(yi, vi, data = dat_sp, method = "REML")
  pred_sp<- predict(fit_sp, transf = transf.ilogit)

  prev   <- sum(data$CT_pos) / sum(data$N)
  Se_p   <- pred_se$pred
  Sp_p   <- pred_sp$pred
  NPV    <- (Sp_p*(1-prev)) / (Sp_p*(1-prev) + (1-Se_p)*prev)

  list(
    data = data, label = label,
    k = nrow(data), N = sum(data$N), CT_pos = sum(data$CT_pos), prev = prev,
    Se = list(est = Se_p, lo = pred_se$ci.lb, hi = pred_se$ci.ub,
              I2 = fit_se$I2, tau2 = fit_se$tau2,
              Q = fit_se$QE, df = fit_se$k-1, pQ = fit_se$QEp,
              fit = fit_se, dat = dat_se),
    Sp = list(est = Sp_p, lo = pred_sp$ci.lb, hi = pred_sp$ci.ub,
              I2 = fit_sp$I2, tau2 = fit_sp$tau2,
              Q = fit_sp$QE, df = fit_sp$k-1, pQ = fit_sp$QEp,
              fit = fit_sp, dat = dat_sp),
    NPV             = NPV,
    CT_avoided      = Sp_p * (1-prev) * 100,
    missed_per_1000 = (1-Se_p) * prev * 1000
  )
}

GFAP_res  <- run_meta_analysis(GFAP_data,  "GFAP+UCH-L1")
S100B_res <- run_meta_analysis(S100B_data, "S100B")

# =============================================================================
# FOREST PLOT helper
# =============================================================================

draw_panel <- function(vals, lo, hi, w, pooled, pooled_lo, pooled_hi,
                       study_names, col_main, xlabel, show_labels, k_label = NULL) {

  k     <- length(vals)
  y_pos <- (k + 0.5):1.5
  pooled_label <- if (!is.null(k_label)) paste0("Pooled [k=", k_label, "]") else "Pooled"

  plot(0, 0, type = "n",
       xlim = c(0, 1.18), ylim = c(0.5, k + 2.0),
       xlab = xlabel, ylab = "", yaxt = "n", xaxt = "n",
       main = "")
  axis(1, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2))

  for (i in 1:k) {
    if (i %% 2 == 1)
      rect(0, i-0.4, 1.18, i+0.4, col = "#F7F7F7", border = NA)
  }
  rect(0, 0.5, 1.18, 1.3,
       col = ifelse(col_main == "#1B4F8A", "#EEF4FF", "#EEF9EE"), border = NA)
  abline(v = 0, col = "#CCCCCC", lwd = 0.5)

  for (i in 1:k) {
    segments(lo[i], y_pos[i], hi[i], y_pos[i], col = col_main, lwd = 1.5)
    segments(lo[i], y_pos[i]-0.1, lo[i], y_pos[i]+0.1, col = col_main, lwd = 1.5)
    segments(hi[i], y_pos[i]-0.1, hi[i], y_pos[i]+0.1, col = col_main, lwd = 1.5)
    pt_size <- 0.5 + (w[i] / max(w)) * 1.5
    points(vals[i], y_pos[i], pch = 21,
           bg = col_main, col = "white", cex = pt_size, lwd = 1)
    text(1.10, y_pos[i], sprintf("%.1f", vals[i]*100),
         adj = 0, cex = 0.75, col = "#555555")
  }

  dw <- (pooled_hi - pooled_lo) / 2
  polygon(c(pooled-dw, pooled, pooled+dw, pooled),
          c(0.9, 1.15, 0.9, 0.65),
          col = col_main, border = "white", lwd = 1.5)
  abline(v = pooled, col = col_main, lwd = 0.8, lty = 2)
  text(1.10, 0.9, sprintf("%.1f", pooled*100),
       adj = 0, cex = 0.85, font = 2, col = col_main)

  if (show_labels) {
    axis(2, at = y_pos,  labels = study_names, las = 1, tick = FALSE, cex.axis = 0.8)
    axis(2, at = 0.9,    labels = pooled_label, las = 1, tick = FALSE, cex.axis = 0.8, font = 2)
  }

  segments(0, 1.35, 1.18, 1.35, col = "#AAAAAA", lwd = 0.8)
}

# =============================================================================
# FIGURE 2 & 3: Forest Plots
# =============================================================================

plot_forest <- function(res, col_main, filename, fig_num) {

  k      <- res$k
  se_n   <- res$data$TP + res$data$FN
  sp_n   <- res$data$TN + res$data$FP
  se_lo  <- sapply(1:k, function(i) wilson_ci(res$data$Se[i], se_n[i])["lo"])
  se_hi  <- sapply(1:k, function(i) wilson_ci(res$data$Se[i], se_n[i])["hi"])
  sp_lo  <- sapply(1:k, function(i) wilson_ci(res$data$Sp[i], sp_n[i])["lo"])
  sp_hi  <- sapply(1:k, function(i) wilson_ci(res$data$Sp[i], sp_n[i])["hi"])
  w_se   <- 1 / (res$data$Se * (1 - res$data$Se) / se_n)
  w_sp   <- 1 / (res$data$Sp * (1 - res$data$Sp) / sp_n)

  N_orig <- sum(res$data$TP - ifelse(res$data$TP %% 1 != 0, 0.5, 0) +
                res$data$FP - ifelse(res$data$FP %% 1 != 0, 0.5, 0) +
                res$data$FN - ifelse(res$data$FN %% 1 != 0, 0.5, 0) +
                res$data$TN - ifelse(res$data$TN %% 1 != 0, 0.5, 0))
  N_orig <- round(N_orig)

  pdf(filename, width = 14, height = max(6, k * 0.75 + 3))
  par(mfrow = c(1, 2), mar = c(5, 12, 2, 2), oma = c(0, 0, 3, 0))

  draw_panel(res$data$Se, se_lo, se_hi, w_se,
             res$Se$est, res$Se$lo, res$Se$hi,
             res$data$study, col_main,
             "Sensitivity (95% CI)", show_labels = TRUE, k_label = res$k)
  text(1.10, 0.3,
       sprintf("I\u00b2=%.0f%%  tau2=%.3f  Q=%.1f (p=%.3f)",
               res$Se$I2, res$Se$tau2, res$Se$Q, res$Se$pQ),
       adj = 1, cex = 0.65, col = "#666666", font = 3)

  draw_panel(res$data$Sp, sp_lo, sp_hi, w_sp,
             res$Sp$est, res$Sp$lo, res$Sp$hi,
             res$data$study, col_main,
             "Specificity (95% CI)", show_labels = FALSE, k_label = NULL)
  text(1.10, 0.3,
       sprintf("I\u00b2=%.0f%%  tau2=%.3f  Q=%.1f (p=%.3f)",
               res$Sp$I2, res$Sp$tau2, res$Sp$Q, res$Sp$pQ),
       adj = 1, cex = 0.65, col = "#666666", font = 3)

  mtext(sprintf("Figure %d. Forest plot: %s  (k=%d, N=%s)",
                fig_num, res$label, res$k, format(N_orig, big.mark=",")),
        side = 3, line = 1.5, outer = TRUE, font = 2, cex = 1.1)
  mtext(sprintf("Pooled: Se = %.1f%% [%.1f\u2013%.1f%%]   Sp = %.1f%% [%.1f\u2013%.1f%%]",
                res$Se$est*100, res$Se$lo*100, res$Se$hi*100,
                res$Sp$est*100, res$Sp$lo*100, res$Sp$hi*100),
        side = 3, line = 0.3, outer = TRUE, cex = 0.85)

  dev.off()
  cat("Saved:", filename, "\n")
}

plot_forest(GFAP_res,  "#1B4F8A", "Fig2_forest_gfap_v15.pdf",  2)
plot_forest(S100B_res, "#1A6B3C", "Fig3_forest_s100b_v15.pdf", 3)

# =============================================================================
# FIGURE 4: SROC Curves
# =============================================================================

plot_sroc <- function(gfap, s100b, filename) {

  pdf(filename, width = 8, height = 8)
  par(mar = c(5, 5, 4, 2))

  plot(0, 0, type = "n", xlim = c(0,1), ylim = c(0,1),
       xlab = "1 \u2212 Specificity (False Positive Rate)",
       ylab = "Sensitivity (True Positive Rate)",
       main = "")

  grid(lty = 1, col = "#EEEEEE")
  abline(0, 1, lty = 3, col = "#AAAAAA")

  max_N <- max(gfap$data$N, s100b$data$N)

  for (res in list(
    list(r = gfap,  col = "#1B4F8A", colL = "#90B8D8", pch = 23,
         fit = suppressWarnings(mada::reitsma(gfap$data[,c("TP","FP","FN","TN")]))),
    list(r = s100b, col = "#1A6B3C", colL = "#90C9A8", pch = 22,
         fit = suppressWarnings(mada::reitsma(s100b$data[,c("TP","FP","FN","TN")])))
  )) {
    sz <- sqrt(res$r$data$N / max_N) * 4
    points(1 - res$r$data$Sp, res$r$data$Se,
           pch = 21, bg = res$colL, col = res$col, cex = sz, lwd = 1.5)
    points(1 - res$r$Sp$est, res$r$Se$est,
           pch = res$pch, bg = res$col, col = "white", cex = 2, lwd = 2)

    mada::sroc(res$fit, add = TRUE, sroclty = 1, sroccol = res$col, lwd = 2.5)
    mada::ROCellipse(res$fit, add = TRUE, col = res$col, lty = 2, lwd = 1)
  }

  legend("bottomright",
    legend = c(
      sprintf("GFAP+UCH-L1  Se=%.1f%%, Sp=%.1f%%",
              gfap$Se$est*100, gfap$Sp$est*100),
      sprintf("S100B          Se=%.1f%%, Sp=%.1f%%",
              s100b$Se$est*100, s100b$Sp$est*100),
      "Chance line (AUC=0.5)",
      "Circle area proportional to N",
      "Dashed ellipse = 95% CI region"
    ),
    col   = c("#1B4F8A","#1A6B3C","#AAAAAA","black","grey40"),
    lty   = c(1, 1, 3, NA, 2),
    lwd   = c(2.5, 2.5, 1, NA, 1),
    pch   = c(23, 22, NA, 21, NA),
    pt.bg = c("#1B4F8A","#1A6B3C", NA,"gray90", NA),
    pt.cex= c(1.5, 1.5, NA, 1, NA),
    bty = "o", bg = "white", box.col = "#CCCCCC", cex = 0.82
  )

  mtext("Figure 4. Summary ROC curves: GFAP+UCH-L1 vs S100B in mTBI",
        side = 1, line = 4, outer = FALSE, font = 2, cex = 0.9)
  mtext(sprintf("GFAP+UCH-L1: Se=%.1f%%, Sp=%.1f%%  |  S100B: Se=%.1f%%, Sp=%.1f%%",
                gfap$Se$est*100, gfap$Sp$est*100,
                s100b$Se$est*100, s100b$Sp$est*100),
        side = 1, line = 5, outer = FALSE, cex = 0.78, col = "#444444")

  dev.off()
  cat("Saved:", filename, "\n")
}

plot_sroc(GFAP_res, S100B_res, "Fig4_sroc_v15.pdf")

# =============================================================================
# FIGURE 5: Clinical Summary Panel
# =============================================================================

plot_clinical_summary <- function(gfap, s100b, filename) {

  col_g  <- "#1B4F8A"; col_s  <- "#1A6B3C"
  col_gl <- "#90B8D8"; col_sl <- "#90C9A8"

  pdf(filename, width = 14, height = 10)
  par(mfrow = c(2,2), mar = c(5,5,3,2), oma = c(0,0,3,0))

  # PANEL A
  se_v  <- c(gfap$Se$est, s100b$Se$est)*100
  sp_v  <- c(gfap$Sp$est, s100b$Sp$est)*100
  se_hi <- c(gfap$Se$hi,  s100b$Se$hi)*100
  se_lo <- c(gfap$Se$lo,  s100b$Se$lo)*100
  sp_hi <- c(gfap$Sp$hi,  s100b$Sp$hi)*100
  sp_lo <- c(gfap$Sp$lo,  s100b$Sp$lo)*100

  bp <- barplot(rbind(se_v, sp_v), beside = TRUE,
                names.arg = c("GFAP+UCH-L1","S100B"),
                col = c(col_g, col_gl, col_s, col_sl),
                ylim = c(0,112), ylab = "Value (%)",
                main = "A.  Sensitivity and Specificity",
                legend.text = c("Sensitivity","Specificity"),
                args.legend = list(x = "topright", cex = 0.85))
  arrows(bp[1,], se_v, bp[1,], se_hi, angle=90, length=0.05, lwd=1.5)
  arrows(bp[1,], se_v, bp[1,], se_lo, angle=90, length=0.05, lwd=1.5)
  arrows(bp[2,], sp_v, bp[2,], sp_hi, angle=90, length=0.05, lwd=1.5)
  arrows(bp[2,], sp_v, bp[2,], sp_lo, angle=90, length=0.05, lwd=1.5)
  text(bp[1,], se_v+3, sprintf("%.1f%%", se_v), cex=0.9, font=2)
  text(bp[2,], sp_v+3, sprintf("%.1f%%", sp_v), cex=0.9, font=2)
  abline(h=95, col="#CC4444", lty=2, lwd=0.8)
  text(max(bp)+0.5, 96, "95%", col="#CC4444", cex=0.75)

  # PANEL B
  npv_v <- c(gfap$NPV, s100b$NPV)*100
  bpb <- barplot(npv_v, names.arg=c("GFAP+UCH-L1","S100B"),
                 col=c(col_g, col_s), ylim=c(97,100.6),
                 ylab="Negative Predictive Value (%)", main="B.  NPV")
  text(bpb, npv_v+0.1, sprintf("%.2f%%", npv_v), cex=1, font=2)
  abline(h=99, col="#CC4444", lty=2, lwd=0.8)
  text(max(bpb)+0.3, 99.1, "99%", col="#CC4444", cex=0.75)

  # PANEL C
  ct_v <- c(gfap$CT_avoided, s100b$CT_avoided)
  ms_v <- c(gfap$missed_per_1000, s100b$missed_per_1000)
  bp2 <- barplot(rbind(ct_v, ms_v), beside=TRUE,
                 names.arg=c("GFAP+UCH-L1","S100B"),
                 col=c(col_g, col_gl, col_s, col_sl),
                 ylab="Value",
                 main="C.  CT Scans Avoided and Missed Lesions",
                 legend.text=c("CT avoided (%)","Missed/1000 pts"),
                 args.legend=list(x="topright", cex=0.85))
  text(bp2[1,], ct_v+0.8, sprintf("%.1f",ct_v), cex=0.9, font=2)
  text(bp2[2,], ms_v+0.1, sprintf("%.1f",ms_v), cex=0.9, font=2)

  # PANEL D
  plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10),
       xlab="", ylab="", axes=FALSE, main="D.  Summary Table")

  rect(0.3, 0.5, 9.7, 9.5, col="white", border="#CCCCCC", lwd=1)
  rect(0.3, 8.5, 9.7, 9.5, col="#E8E8E8", border="#CCCCCC")
  text(c(2.5,5.5,8.2), 9,
       c("Parameter","GFAP+UCH-L1","S100B"),
       font=2, cex=0.85, col=c("black",col_g,col_s))

  rows <- list(
    c("Studies (k)",      as.character(gfap$k), as.character(s100b$k)),
    c("Patients (N)",     format(round(gfap$N,0), big.mark=","), format(round(s100b$N,0), big.mark=",")),
    c("CT+ events",       paste0(round(gfap$CT_pos,0)," (",round(gfap$prev*100,1),"%)"),
                           paste0(round(s100b$CT_pos,0)," (",round(s100b$prev*100,1),"%)")),
    c("Sensitivity (%)",
      paste0(round(gfap$Se$est*100,1)," [",round(gfap$Se$lo*100,1),"\u2013",round(gfap$Se$hi*100,1),"]"),
      paste0(round(s100b$Se$est*100,1)," [",round(s100b$Se$lo*100,1),"\u2013",round(s100b$Se$hi*100,1),"]")),
    c("Specificity (%)",
      paste0(round(gfap$Sp$est*100,1)," [",round(gfap$Sp$lo*100,1),"\u2013",round(gfap$Sp$hi*100,1),"]"),
      paste0(round(s100b$Sp$est*100,1)," [",round(s100b$Sp$lo*100,1),"\u2013",round(s100b$Sp$hi*100,1),"]")),
    c("I\u00b2 Se / Sp (%)",
      paste0(round(gfap$Se$I2,0)," / ",round(gfap$Sp$I2,0)),
      paste0(round(s100b$Se$I2,0)," / ",round(s100b$Sp$I2,0))),
    c("NPV (%)",          round(gfap$NPV*100,2),  round(s100b$NPV*100,2)),
    c("CT avoided (%)",   round(gfap$CT_avoided,1),round(s100b$CT_avoided,1)),
    c("Missed/1000 pts",  round(gfap$missed_per_1000,1),round(s100b$missed_per_1000,1))
  )

  yp <- seq(7.8, 1.0, length.out=length(rows))
  for (i in seq_along(rows)) {
    if (i%%2==0) rect(0.3,yp[i]-0.38,9.7,yp[i]+0.38,col="#F8F8F8",border=NA)
    text(c(2.5,5.5,8.2), yp[i], rows[[i]], cex=0.75)
    segments(0.3,yp[i]-0.44,9.7,yp[i]-0.44,col="#DDDDDD",lwd=0.5)
  }

  mtext("Figure 5.  Clinical Summary \u2014 GFAP+UCH-L1 vs S100B for mTBI Triage",
        side=3, line=1.5, outer=TRUE, font=2, cex=1.1)
  mtext("Dataset v15 definitivo | Bivariate random-effects (REML) | metafor + mada",
        side=3, line=0.3, outer=TRUE, cex=0.85)

  dev.off()
  cat("Saved:", filename, "\n")
}

plot_clinical_summary(GFAP_res, S100B_res, "Fig5_clinical_summary_v15.pdf")

# =============================================================================
# FIGURE S1: Funnel Plots
# =============================================================================

plot_funnel <- function(res, col_main, filename) {

  pdf(filename, width=10, height=6)
  par(mfrow=c(1,2), mar=c(5,4,3,2), oma=c(0,0,3,0))

  n_studies <- length(res$data$study)
  pos_alt <- rep(3, n_studies)
  pos_alt[seq(2, n_studies, by = 2)] <- 1

  funnel(res$Se$fit,
         main     = "Sensitivity",
         xlab     = "Logit(Sensitivity)",
         col      = col_main,
         bg       = col_main)
  text(res$Se$fit$yi, sqrt(res$Se$fit$vi), labels = res$data$study,
       cex = 0.5, pos = pos_alt, offset = 0.35, col = "grey30")

  funnel(res$Sp$fit,
         main     = "Specificity",
         xlab     = "Logit(Specificity)",
         col      = col_main,
         bg       = col_main)
  text(res$Sp$fit$yi, sqrt(res$Sp$fit$vi), labels = res$data$study,
       cex = 0.5, pos = pos_alt, offset = 0.35, col = "grey30")

  mtext(paste("Figure S1.  Funnel Plots \u2014", res$label),
        side=3, line=1.2, outer=TRUE, font=2, cex=1.05)
  mtext("Egger test p-values: see Results section",
        side=3, line=0.1, outer=TRUE, cex=0.82)

  dev.off()
  cat("Saved:", filename, "\n")
}

plot_funnel(GFAP_res,  "#1B4F8A", "FigS1_funnel_gfap_v15.pdf")
plot_funnel(S100B_res, "#1A6B3C", "FigS1_funnel_s100b_v15.pdf")

# =============================================================================
# EGGER TEST
# =============================================================================

cat("\n", rep("=", 72), "\n", sep="")
cat(" PUBLICATION BIAS \u2014 Egger Test\n")
cat(rep("=", 72), "\n", sep="")

for (res in list(GFAP_res, S100B_res)) {
  prec     <- 1 / sqrt(res$Se$dat$vi)
  std_eff  <- (res$Se$dat$yi - log(res$Se$est/(1-res$Se$est))) / sqrt(res$Se$dat$vi)
  egger    <- lm(std_eff ~ prec)
  pval     <- summary(egger)$coefficients[1, 4]
  cat(sprintf("\n%s (Sensitivity):\n  Intercept = %.3f  |  p = %.3f  |  %s\n",
              res$label, coef(egger)[1], pval,
              ifelse(pval > 0.05, "No evidence of publication bias",
                     "Possible publication bias")))
}

cat("\n NOTE: Egger test ha bassa potenza con k<10. Valutare assieme al funnel plot.\n")
cat(rep("=", 72), "\n", sep="")

# =============================================================================
# CONSOLE OUTPUT
# =============================================================================

print_summary <- function(res) {
  cat("\n", rep("=",72), "\n", sep="")
  cat(" ", res$label, "\u2014 Bivariate RE (metafor REML)\n")
  cat(rep("=",72), "\n", sep="")
  cat(sprintf(" k=%d  N=%d  CT+=%d (%.1f%%)\n",
              res$k, round(res$N,0), round(res$CT_pos,0), res$prev*100))
  cat(rep("-",72), "\n", sep="")

  cat(sprintf(" %-30s %6s %16s %6s %16s\n",
              "Study","Se%","[95%CI]","Sp%","[95%CI]"))
  cat(rep("-",72), "\n", sep="")

  se_n <- res$data$TP + res$data$FN
  sp_n <- res$data$TN + res$data$FP
  for (i in 1:res$k) {
    ci_se <- wilson_ci(res$data$Se[i], se_n[i])
    ci_sp <- wilson_ci(res$data$Sp[i], sp_n[i])
    cat(sprintf(" %-30s %6.1f [%5.1f-%5.1f] %6.1f [%5.1f-%5.1f]\n",
                res$data$study[i],
                res$data$Se[i]*100, ci_se["lo"]*100, ci_se["hi"]*100,
                res$data$Sp[i]*100, ci_sp["lo"]*100, ci_sp["hi"]*100))
  }
  cat(rep("-",72), "\n", sep="")

  cat(sprintf("\n POOLED (REML):\n"))
  cat(sprintf("   Sensitivity  %.1f%% [%.1f%% \u2013 %.1f%%]\n",
              res$Se$est*100, res$Se$lo*100, res$Se$hi*100))
  cat(sprintf("   Specificity  %.1f%% [%.1f%% \u2013 %.1f%%]\n",
              res$Sp$est*100, res$Sp$lo*100, res$Sp$hi*100))

  cat(sprintf("\n HETEROGENEITY:\n"))
  cat(sprintf("   Se: I\u00b2=%.0f%%  \u03c4\u00b2=%.4f  Q=%.2f (df=%d, p=%.4f)\n",
              res$Se$I2, res$Se$tau2, res$Se$Q, res$Se$df, res$Se$pQ))
  cat(sprintf("   Sp: I\u00b2=%.0f%%  \u03c4\u00b2=%.4f  Q=%.2f (df=%d, p=%.4f)\n",
              res$Sp$I2, res$Sp$tau2, res$Sp$Q, res$Sp$df, res$Sp$pQ))

  cat(sprintf("\n CLINICAL METRICS (prev=%.1f%%):\n", res$prev*100))
  cat(sprintf("   NPV                  %.2f%%\n", res$NPV*100))
  cat(sprintf("   CT scans avoided     %.1f%%\n", res$CT_avoided))
  cat(sprintf("   Missed lesions/1000  %.1f\n",   res$missed_per_1000))
  cat(rep("=",72), "\n", sep="")
}

print_summary(GFAP_res)
print_summary(S100B_res)

cat("\n Files generated:\n")
for (f in c("Fig2_forest_gfap_v15.pdf","Fig3_forest_s100b_v15.pdf",
            "Fig4_sroc_v15.pdf","Fig5_clinical_summary_v15.pdf",
            "FigS1_funnel_gfap_v15.pdf","FigS1_funnel_s100b_v15.pdf"))
  cat("  ", f, "\n")
