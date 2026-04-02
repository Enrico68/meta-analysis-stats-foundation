# ============================================================
# mTBI Biomarkers DTA Meta-Analysis — v7
# ============================================================

if (!require("mada")) install.packages("mada", repos="https://cran.r-project.org")
if (!require("metafor")) install.packages("metafor", repos="https://cran.r-project.org")
library(mada); library(metafor)

# === DATA ===
gfap <- data.frame(
  study=c("Bazarian 2021","Lagares 2024","Faisal 2023","Legramante 2024",
          "Ladang 2025","Milevoj K. 2025","Curran 2025","Campagna 2025"),
  TP=c(115,176,12,7,112,107,12,8), FP=c(1061,946,75,89,148,529,56,59),
  FN=c(5,3,0,0,1,5,0,0), TN=c(720,313,43,34,101,181,21,30))
gfap$N <- gfap$TP+gfap$FP+gfap$FN+gfap$TN

s100b <- data.frame(
  study=c("Oris 2024","Oris 2021","Hopman 2023"),
  TP=c(7,63,69), FP=c(148,891,355), FN=c(0,0,5), TN=c(51,218,66))
s100b$N <- s100b$TP+s100b$FP+s100b$FN+s100b$TN

gfap_sens <- data.frame(
  study=c("Lapic 2025","Jalali 2025","Spaziani 2025","Welch 2025"),
  TP=c(27,5,72,116), FP=c(104,83,277,1066), FN=c(2,1,3,4), TN=c(42,34,89,713))
gfap_sens$N <- gfap_sens$TP+gfap_sens$FP+gfap_sens$FN+gfap_sens$TN

# === REITSMA MODELS ===
cat("\n========== GFAP+UCH-L1 (8 studies, N=",sum(gfap$N),") ==========\n\n")
fit_gfap <- reitsma(gfap, correction.control="all", correction=0.5)
print(summary(fit_gfap))
npv <- sum(gfap$TN)/(sum(gfap$TN)+sum(gfap$FN))
ct_red <- (sum(gfap$FN)+sum(gfap$TN))/sum(gfap$N)*100
cat(sprintf("\nNPV=%.3f | CT-red=%.1f%% | Missed=%d/%d (%.1f%%)\n",
            npv,ct_red,sum(gfap$FN),sum(gfap$TP+gfap$FN),
            sum(gfap$FN)/sum(gfap$TP+gfap$FN)*100))

cat("\n\n========== S100B (3 studies, N=",sum(s100b$N),") ==========\n\n")
s100_ok <- FALSE
tryCatch({
  fit_s100 <- reitsma(s100b, correction.control="all", correction=0.5)
  print(summary(fit_s100)); s100_ok <- TRUE
}, error=function(e) cat("Did not converge\n"))
npv_s <- sum(s100b$TN)/(sum(s100b$TN)+sum(s100b$FN))
ct_red_s <- (sum(s100b$FN)+sum(s100b$TN))/sum(s100b$N)*100
cat(sprintf("\nNPV=%.3f | CT-red=%.1f%% | Missed=%d/%d (%.1f%%)\n",
            npv_s,ct_red_s,sum(s100b$FN),sum(s100b$TP+s100b$FN),
            sum(s100b$FN)/sum(s100b$TP+s100b$FN)*100))

# === FOREST PLOTS (manual, no rect) ===
cat("\n--- Generating Plots ---\n")

make_forest <- function(data, fit, title, filename, col_main="steelblue") {
  cc <- ifelse(data$FN==0 | data$TP==0, 0.5, 0)
  tp <- data$TP+cc; fn <- data$FN+cc; fp <- data$FP+cc; tn <- data$TN+cc
  se <- tp/(tp+fn); sp <- tn/(tn+fp)
  k <- nrow(data)
  
  wilson_ci <- function(x, n) {
    if (n==0) return(c(0,1))
    z <- 1.96; p <- x/n; d <- 1+z^2/n
    mid <- p+z^2/(2*n); off <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))
    c(max(0,(mid-off)/d), min(1,(mid+off)/d))
  }
  se_ci <- t(sapply(1:k, function(i) wilson_ci(data$TP[i], data$TP[i]+data$FN[i])))
  sp_ci <- t(sapply(1:k, function(i) wilson_ci(data$TN[i], data$TN[i]+data$FP[i])))
  
  pooled_se <- fit$Sens; se_pool_ci <- fit$Sens.ci
  pooled_sp <- 1-fit$FPR; sp_pool_ci <- c(1-fit$FPR.ci[2], 1-fit$FPR.ci[1])
  
  png(filename, width=1400, height=250+k*80, res=150)
  par(mfrow=c(1,2), mar=c(5,11,4,3))
  y <- k:1
  
  # Sensitivity
  plot(NA, xlim=c(0.75,1.02), ylim=c(0,k+1.5), xlab="Sensitivity", ylab="", yaxt="n",
       main=paste0("Sensitivity\nPooled: ",sprintf("%.1f%% (%.1f-%.1f)",
                   pooled_se*100,se_pool_ci[1]*100,se_pool_ci[2]*100)))
  abline(v=pooled_se, col="red", lty=2, lwd=1.5)
  for (i in 1:k) {
    segments(se_ci[i,1],y[i],se_ci[i,2],y[i], col=col_main, lwd=2)
    points(se[i],y[i], pch=15, cex=1.5, col=col_main)
    text(1.018, y[i], sprintf("%.1f%%",se[i]*100), cex=0.65, xpd=TRUE)
  }
  points(pooled_se, 0.3, pch=18, cex=2.5, col="red")
  segments(se_pool_ci[1],0.3,se_pool_ci[2],0.3, col="red", lwd=2)
  axis(2, at=y, labels=data$study, las=1, cex.axis=0.8)
  axis(2, at=0.3, labels="Pooled", las=1, cex.axis=0.8, font=2, col.axis="red")
  
  # Specificity
  par(mar=c(5,2,4,9))
  plot(NA, xlim=c(0,0.55), ylim=c(0,k+1.5), xlab="Specificity", ylab="", yaxt="n",
       main=paste0("Specificity\nPooled: ",sprintf("%.1f%% (%.1f-%.1f)",
                   pooled_sp*100,sp_pool_ci[1]*100,sp_pool_ci[2]*100)))
  abline(v=pooled_sp, col="red", lty=2, lwd=1.5)
  for (i in 1:k) {
    segments(sp_ci[i,1],y[i],sp_ci[i,2],y[i], col=col_main, lwd=2)
    points(sp[i],y[i], pch=15, cex=1.5, col=col_main)
    text(0.545, y[i], sprintf("%.1f%%",sp[i]*100), cex=0.65, xpd=TRUE)
  }
  points(pooled_sp, 0.3, pch=18, cex=2.5, col="red")
  segments(sp_pool_ci[1],0.3,sp_pool_ci[2],0.3, col="red", lwd=2)
  
  mtext(title, side=3, outer=FALSE, line=-1, cex=0.8, font=2, at=-0.15)
  dev.off()
}

make_forest(gfap, fit_gfap, "GFAP+UCH-L1 (Reitsma Bivariate)",
            "forest_GFAP_UCH_L1_reitsma.png", "steelblue")
cat("Saved: forest_GFAP_UCH_L1_reitsma.png\n")

if (s100_ok) {
  make_forest(s100b, fit_s100, "S100B (Reitsma Bivariate)",
              "forest_S100B_reitsma.png", "forestgreen")
  cat("Saved: forest_S100B_reitsma.png\n")
}

# === SROC ===
png("SROC_GFAP_reitsma.png", width=900, height=900, res=150)
plot(fit_gfap, sroclwd=2, predict=TRUE, conf.region=TRUE,
     main="SROC - GFAP+UCH-L1 (Bivariate Reitsma)")
fpr_g <- 1-gfap$TN/(gfap$TN+gfap$FP)
tpr_g <- gfap$TP/(gfap$TP+gfap$FN+0.5*(gfap$FN==0))
points(fpr_g,tpr_g, pch=19, cex=sqrt(gfap$N)/25, col=rgb(0.1,0.3,0.8,0.6))
text(fpr_g,tpr_g, labels=gsub(" \\d+$","",gfap$study), pos=4, cex=0.55, col="darkblue")
dev.off()
cat("Saved: SROC_GFAP_reitsma.png\n")

png("SROC_comparison_reitsma.png", width=1000, height=1000, res=150)
plot(fit_gfap, sroclwd=2, col="blue", predict=FALSE, conf.region=TRUE,
     main="SROC: GFAP+UCH-L1 vs S100B", xlim=c(0,1), ylim=c(0.8,1.02))
points(fpr_g,tpr_g, pch=19, cex=sqrt(gfap$N)/25, col=rgb(0.1,0.3,0.8,0.6))
fpr_s <- 1-s100b$TN/(s100b$TN+s100b$FP)
tpr_s <- s100b$TP/(s100b$TP+s100b$FN+0.5*(s100b$FN==0))
points(fpr_s,tpr_s, pch=17, cex=sqrt(s100b$N)/25, col=rgb(0.1,0.5,0.1,0.7))
text(fpr_s,tpr_s, labels=s100b$study, pos=4, cex=0.55, col="darkgreen")
legend("bottomright", legend=c("GFAP+UCH-L1","S100B"),
       col=c("blue","darkgreen"), pch=c(19,17), cex=0.8, bty="n")
dev.off()
cat("Saved: SROC_comparison_reitsma.png\n")

# === DEEKS ===
cat("\n========== DEEKS FUNNEL ==========\n\n")
gfap$DOR <- ((gfap$TP+0.5)*(gfap$TN+0.5))/((gfap$FP+0.5)*(gfap$FN+0.5))
gfap$lnDOR <- log(gfap$DOR)
gfap$var_lnDOR <- 1/(gfap$TP+0.5)+1/(gfap$FP+0.5)+1/(gfap$FN+0.5)+1/(gfap$TN+0.5)
gfap$ESS <- 4*gfap$N/(((gfap$TP+gfap$FN)/gfap$N)+((gfap$FP+gfap$TN)/gfap$N))^2
gfap$inv_sqrt_ESS <- 1/sqrt(gfap$ESS)
deeks <- lm(lnDOR~inv_sqrt_ESS, data=gfap, weights=1/var_lnDOR)
deeks_p <- summary(deeks)$coefficients["inv_sqrt_ESS","Pr(>|t|)"]
cat(sprintf("Deeks p=%.3f",deeks_p))
if(deeks_p>0.10) cat(" => No bias\n") else cat(" => Possible bias\n")

png("deeks_funnel_GFAP.png", width=900, height=700, res=150)
plot(gfap$inv_sqrt_ESS, gfap$lnDOR, pch=19, cex=sqrt(gfap$N)/25, col="steelblue",
     xlab="1/sqrt(ESS)", ylab="ln(DOR)",
     main=sprintf("Deeks Funnel - GFAP+UCH-L1 (p=%.3f)",deeks_p))
abline(deeks, col="red", lwd=2, lty=2)
text(gfap$inv_sqrt_ESS, gfap$lnDOR, labels=gfap$study, pos=4, cex=0.6); grid()
dev.off()
cat("Saved: deeks_funnel_GFAP.png\n")

# === SENSITIVITY ANALYSES ===
cat("\n========== SENSITIVITY ANALYSES ==========\n\n")
run_sa <- function(data,label) {
  tryCatch({
    fit <- reitsma(data, correction.control="all", correction=0.5)
    cat(sprintf("%s: Se=%.3f(%.3f-%.3f) Sp=%.3f(%.3f-%.3f) AUC=%.3f [k=%d N=%d]\n",
                label,fit$Sens,fit$Sens.ci[1],fit$Sens.ci[2],
                fit$Spec,fit$Spec.ci[1],fit$Spec.ci[2],
                summary(fit)$AUC,nrow(data),sum(data$N)))
  }, error=function(e) cat(sprintf("%s: Did not converge\n",label)))
}
run_sa(rbind(gfap[,c("TP","FP","FN","TN","N")],
             gfap_sens[gfap_sens$study=="Lapic 2025",c("TP","FP","FN","TN","N")]),"SA1 +Lapic")
run_sa(rbind(gfap[,c("TP","FP","FN","TN","N")],
             gfap_sens[gfap_sens$study!="Welch 2025",c("TP","FP","FN","TN","N")]),"SA2 All-Welch")
run_sa(rbind(gfap[gfap$study!="Bazarian 2021",c("TP","FP","FN","TN","N")],
             data.frame(TP=116,FP=1066,FN=4,TN=713,N=1899)),"SA3 Welch/Baz")
run_sa(rbind(s100b[,c("TP","FP","FN","TN","N")],
             data.frame(TP=11,FP=108,FN=0,TN=60,N=179)),"SA4 S100B+Kah")

# === FINAL ===
cat("\n========== FINAL RESULTS ==========\n\n")
cat(sprintf("GFAP+UCH-L1: Se=%.1f%%(%.1f-%.1f) Sp=%.1f%%(%.1f-%.1f) AUC=%.3f NPV=%.1f%% CT-red=%.1f%%\n",
            fit_gfap$Sens*100,fit_gfap$Sens.ci[1]*100,fit_gfap$Sens.ci[2]*100,
            (1-fit_gfap$FPR)*100,(1-fit_gfap$FPR.ci[2])*100,(1-fit_gfap$FPR.ci[1])*100,
            summary(fit_gfap)$AUC,npv*100,ct_red))
if(s100_ok) cat(sprintf("S100B:       Se=%.1f%%(%.1f-%.1f) Sp=%.1f%%(%.1f-%.1f) pAUC=%.3f NPV=%.1f%% CT-red=%.1f%%\n",
            fit_s100$Sens*100,fit_s100$Sens.ci[1]*100,fit_s100$Sens.ci[2]*100,
            (1-fit_s100$FPR)*100,(1-fit_s100$FPR.ci[2])*100,(1-fit_s100$FPR.ci[1])*100,
            summary(fit_s100)$pAUC,npv_s*100,ct_red_s))
cat("\n=== COMPLETE ===\n")