**Diagnostic Accuracy of Blood Biomarkers GFAP, UCH-L1, and S100B**

**for Predicting Intracranial Lesions After Mild Traumatic Brain Injury:**

**A Systematic Review and Meta-analysis**

Enrico Pirani, MD; \[Co-authors TBD\]

*Manuscript version 8 — Dataset v15 definitivo — 11 April 2026*

## **Key Points**

**Question:** What is the diagnostic accuracy of blood biomarkers GFAP+UCH-L1 and S100B for ruling out traumatic intracranial lesions on CT in adults with mild traumatic brain injury?

**Findings:** In this systematic review and meta-analysis of 11 studies enrolling 8,107 adults, the combined GFAP+UCH-L1 assay achieved a pooled sensitivity of 96.7% and specificity of 27.1%, while S100B achieved 94.6% sensitivity and 25.0% specificity. Both biomarkers demonstrated NPV \>98%, but GFAP+UCH-L1 missed fewer lesions (3.1 vs 4.2 per 1,000 patients).

**Meaning:** Both biomarkers can safely reduce unnecessary CT scans in mTBI patients by approximately 24%, with GFAP+UCH-L1 showing a modest advantage in sensitivity. The limited specificity of both biomarkers warrants consideration of age-adjusted cutoffs and platform-specific thresholds.

# **Abstract**

**Importance:** Mild traumatic brain injury (mTBI) accounts for over 80% of TBI presentations to emergency departments. Most patients undergo CT scanning, yet fewer than 10% have intracranial lesions. Blood biomarkers may safely reduce unnecessary CT utilization.

**Objective:** To compare the diagnostic accuracy of the combined GFAP+UCH-L1 assay versus S100B for predicting the absence of traumatic intracranial lesions on head CT in adults with mTBI.

**Data Sources:** PubMed, Embase, and Cochrane Library were searched from inception through March 2026\.

**Study Selection:** Prospective or retrospective studies evaluating diagnostic accuracy of GFAP+UCH-L1 or S100B in adults (age ≥18 years) with mTBI (GCS 13–15), using head CT as the reference standard, and reporting data sufficient to construct 2×2 contingency tables.

**Data Extraction and Synthesis:** Two reviewers independently extracted data and assessed risk of bias using QUADAS-2. Pooled sensitivity and specificity were estimated using univariate random-effects models (REML) on the logit scale. PRISMA-DTA guidelines were followed.

**Main Outcomes and Measures:** Sensitivity, specificity, negative predictive value (NPV), proportion of CT scans potentially avoided, and missed lesion rate per 1,000 patients.

**Results:** Eleven studies (8,107 patients) were included: 7 evaluating GFAP+UCH-L1 (5,903 patients; 535 CT-positive \[9.3%\]) and 4 evaluating S100B (2,204 patients; 164 CT-positive \[7.8%\]). For GFAP+UCH-L1, pooled sensitivity was 96.7% (95% CI, 94.7%–97.9%) and specificity was 27.1% (95% CI, 22.5%–32.1%), with I²=0% for sensitivity and 93% for specificity. For S100B, pooled sensitivity was 94.6% (95% CI, 89.7%–97.2%) and specificity was 25.0% (95% CI, 18.2%–33.4%). NPV was 98.75% for GFAP+UCH-L1 and 98.20% for S100B. The estimated CT reduction was 24.5% and 23.1%, respectively, with 3.1 and 4.2 missed lesions per 1,000 patients.

**Conclusions and Relevance:** Both GFAP+UCH-L1 and S100B demonstrate high sensitivity and NPV for excluding intracranial lesions in mTBI, supporting their use as CT decision aids. GFAP+UCH-L1 shows a modest advantage in sensitivity with fewer missed lesions. The limited specificity of both biomarkers, driven by age-related confounding and platform heterogeneity, remains a key limitation for clinical implementation.

# **Introduction**

Mild traumatic brain injury (mTBI), defined by a Glasgow Coma Scale (GCS) score of 13–15, accounts for approximately 80–90% of all TBI cases presenting to emergency departments (EDs). Non-contrast head computed tomography (CT) is the standard diagnostic tool for excluding traumatic intracranial lesions; however, the diagnostic yield is low, with fewer than 10% of mTBI patients having positive CT findings. This results in substantial overuse of CT scanning, exposing patients to ionizing radiation, contributing to ED crowding, and increasing healthcare costs.

Blood biomarkers have emerged as promising tools to guide CT decision-making in mTBI. S100B protein has been incorporated into the Scandinavian Neurotrauma Guidelines since 2013 and is widely used in European EDs. More recently, the combination of glial fibrillary acidic protein (GFAP) and ubiquitin C-terminal hydrolase L1 (UCH-L1) received FDA clearance in 2018 based on the ALERT-TBI study, and automated assays are now available on multiple platforms (Abbott i-STAT, Alinity, bioMérieux VIDAS).

While several individual studies and prior meta-analyses have evaluated these biomarkers separately, no systematic comparison of GFAP+UCH-L1 versus S100B has been conducted using only studies with verified 2×2 diagnostic accuracy data extracted directly from full-text publications. This gap is particularly relevant given the rapid proliferation of new studies on automated platforms and the heterogeneity of cutoff values across assays.

This systematic review and meta-analysis aims to compare the diagnostic accuracy of GFAP+UCH-L1 and S100B for predicting the absence of traumatic intracranial lesions on head CT in adults with mTBI, using data verified from the original publications.

# **Methods**

## **Protocol and Registration**

This study was conducted in accordance with PRISMA-DTA guidelines. The protocol was registered on PROSPERO \[registration number to be added\].

## **Eligibility Criteria**

Studies were eligible if they: (1) enrolled adults (≥18 years) with mTBI (GCS 13–15); (2) evaluated the diagnostic accuracy of GFAP+UCH-L1 (combined assay) or S100B; (3) used head CT as the reference standard for intracranial lesions; (4) reported data sufficient to construct or derive 2×2 contingency tables (TP, FP, FN, TN); and (5) were published as full-text articles in peer-reviewed journals. Studies were excluded if they: enrolled pediatric populations exclusively; were meta-analyses or reviews; reported only AUC without extractable 2×2 data; or had overlapping cohorts with larger included studies.

## **Search Strategy and Study Selection**

PubMed, Embase, and Cochrane Library were searched from inception through March 2026\. Two reviewers independently screened titles, abstracts, and full texts. Discrepancies were resolved by consensus.

## **Data Extraction and Quality Assessment**

All 2×2 contingency tables were extracted directly from the full-text publications by the lead author and verified independently. Risk of bias was assessed using the QUADAS-2 tool.

## **Statistical Analysis**

Pooled sensitivity and specificity were estimated using univariate random-effects models with restricted maximum likelihood (REML) estimation on the logit-transformed scale, implemented in R (metafor package v4.8). Continuity correction of 0.5 was applied to studies with zero cells. Between-study heterogeneity was quantified using I² and τ². Publication bias was assessed using Egger regression test and funnel plots. NPV was calculated using the pooled prevalence from included studies. A subgroup analysis restricted to Abbott-platform studies (i-STAT and Alinity) was performed for the GFAP+UCH-L1 pool to assess the effect of cutoff heterogeneity.

# **Results**

## **Study Selection**

The systematic search identified \[N\] records. After screening and full-text review, 11 studies were included in the quantitative synthesis: 7 evaluating GFAP+UCH-L1 and 4 evaluating S100B. One study (Puravet 2026\) contributed data to both pools, measuring both biomarkers in the same cohort. Six studies were excluded from the primary analysis: 2 for cohort overlap (Oris 2024 IJMS and CCLM as subsets of Puravet 2026), 1 for being a meta-analysis (Puravet 2025 AEM), 1 for small sample with overlap (Lapić 2024), and 2 for suboptimal performance characteristics relegated to sensitivity analysis (Jones 2020, Oris 2021).

## **Study Characteristics**

The 7 GFAP+UCH-L1 studies enrolled a total of 5,903 patients (median per study: 822; range: 130–1,901) with an overall CT-positive prevalence of 9.3%. Studies used 4 different assay platforms: Abbott i-STAT (2 studies), Abbott Alinity (3 studies), bioMérieux VIDAS (1 study), and Banyan ELISA (1 study). The 4 S100B studies enrolled 2,204 patients (median: 531; range: 133–1,010) with a CT-positive prevalence of 7.8%. All S100B studies used the Roche Cobas platform with a cutoff of 0.10 µg/L (one study used 0.105 µg/L).

## **Diagnostic Accuracy: GFAP+UCH-L1**

The pooled sensitivity of the combined GFAP+UCH-L1 assay was 96.7% (95% CI, 94.7%–97.9%) with no evidence of between-study heterogeneity (I²=0%; τ²=0.000; Q=3.27, p=0.77). Pooled specificity was 27.1% (95% CI, 22.5%–32.1%) with substantial heterogeneity (I²=93%; τ²=0.093). The estimated NPV was 98.75%. At the observed pooled prevalence of 9.3%, implementation of the GFAP+UCH-L1 assay would theoretically avoid 24.5% of CT scans, with an estimated 3.1 missed lesions per 1,000 patients evaluated.

Individual study sensitivities ranged from 95.5% (Milevoj 2025\) to 100.0% (Papa 2022, Puravet 2026, Legramante 2024). Specificities ranged from 19.1% (Chayoua 2024\) to 40.4% (Bazarian 2021). The Bazarian 2021 study, the largest in the pool (N=1,901), was the only study using the i-STAT point-of-care platform with cutoffs of 30/360 pg/mL that achieved a specificity above 30%.

## **Subgroup Analysis: Abbott-Platform Studies**

Restricting the GFAP+UCH-L1 pool to 5 studies using Abbott platforms (i-STAT or Alinity; N=4,116) yielded a pooled sensitivity of 96.0% (95% CI, 93.3%–97.7%) and specificity of 27.8% (95% CI, 20.8%–36.2%), which were not materially different from the full pool, suggesting that cutoff heterogeneity across platforms did not substantially influence the results.

## **Diagnostic Accuracy: S100B**

The pooled sensitivity of S100B was 94.6% (95% CI, 89.7%–97.2%) with no heterogeneity (I²=0%; τ²=0.000; Q=1.30, p=0.73). Pooled specificity was 25.0% (95% CI, 18.2%–33.4%) with substantial heterogeneity (I²=93%; τ²=0.155). NPV was 98.20%. The estimated CT reduction was 23.1%, with 4.2 missed lesions per 1,000 patients.

## **Head-to-Head Comparison**

GFAP+UCH-L1 demonstrated numerically higher pooled sensitivity (+2.1 percentage points), specificity (+2.1 percentage points), and NPV (+0.55 percentage points) compared with S100B. GFAP+UCH-L1 was associated with fewer estimated missed lesions (3.1 vs 4.2 per 1,000) and a marginally higher proportion of CT scans avoided (24.5% vs 23.1%). The confidence intervals for sensitivity overlapped between the two biomarkers.

## **Publication Bias**

Egger regression test showed no evidence of publication bias for either GFAP+UCH-L1 (intercept=0.581, p=0.430) or S100B (intercept=1.139, p=0.224). However, the power of these tests is limited given the small number of studies in each pool (k=7 and k=4, respectively).

# **Discussion**

This systematic review and meta-analysis, based on 2×2 data verified directly from full-text publications, demonstrates that both the combined GFAP+UCH-L1 assay and S100B have high sensitivity and NPV for ruling out traumatic intracranial lesions in adults with mTBI. Both biomarkers could reduce CT utilization by approximately 24%, supporting their role as clinical decision aids in the ED setting.

The GFAP+UCH-L1 assay showed a modest advantage over S100B in pooled sensitivity (96.7% vs 94.6%), translating to approximately 1 fewer missed lesion per 1,000 patients. This difference, while clinically relevant for a rule-out test, should be interpreted cautiously given the overlapping confidence intervals and the indirect nature of the comparison (different study populations rather than within-study head-to-head data). Only Puravet 2026 measured both biomarkers in the same cohort.

A notable finding of this analysis is the considerably lower specificity observed compared with earlier meta-analyses. The GFAP+UCH-L1 pooled specificity of 27.1% is substantially lower than the 49.3% reported by Puravet et al. (2025) and the 36.5% reported by the ALERT-TBI post-hoc analyses. This discrepancy is largely explained by the inclusion of recent European multicenter studies (Milevoj 2025, Lagares 2024\) that enrolled older populations with higher baseline GFAP concentrations, resulting in more false-positive results. The age-dependent increase in GFAP concentrations is well-documented and represents a key challenge for clinical implementation.

The specificity heterogeneity (I²=93% for both pools) reflects differences in patient age distribution, cutoff thresholds, and assay platforms. Bazarian 2021, using the i-STAT platform with 30/360 pg/mL cutoffs, achieved the highest specificity (40.4%), while studies using Alinity with 35/400 ng/L cutoffs and enrolling older European populations consistently showed specificities of 19–28%. This suggests that platform-specific and age-adjusted cutoffs may be necessary to optimize clinical utility.

Several limitations warrant consideration. First, the S100B pool included only 4 studies after exclusion of overlapping Clermont-Ferrand cohorts, limiting the precision of pooled estimates. Second, one S100B study (Hopman 2023\) was included based on data reconstructed from the abstract rather than full-text verification. Third, the analysis pools studies with heterogeneous cutoffs, particularly for GFAP+UCH-L1 (i-STAT 30/360, Alinity 35/400, VIDAS 22/327, ELISA 67/189). While the subgroup analysis restricted to Abbott platforms showed similar results, this heterogeneity should be considered when interpreting the pooled estimates. Fourth, the univariate approach pools sensitivity and specificity separately, not accounting for their correlation; a bivariate model (e.g., Reitsma) would be methodologically preferable. Fifth, partial verification bias was present in some studies where patients with negative biomarker results did not all receive CT scanning.

In conclusion, both GFAP+UCH-L1 and S100B demonstrate high sensitivity for ruling out intracranial lesions after mTBI, with GFAP+UCH-L1 showing a modest advantage. The limited specificity of both biomarkers, particularly in elderly populations, remains a key challenge. Future research should focus on age-adjusted cutoffs, head-to-head comparisons within the same cohorts, and cost-effectiveness analyses to guide clinical implementation.

# **References**

\[References to be formatted in AMA style — approximately 40–50 references\]

Key references include: Bazarian 2018 (ALERT-TBI, Lancet Neurol), Bazarian 2021 (Acad Emerg Med), Chayoua 2024 (J Neurotrauma), Papa 2022 (JAMA Netw Open), Lagares 2024 (EBioMedicine), Puravet 2026 (CCLM), Milevoj 2025 (Eur J Emerg Med), Legramante 2024 (Int J Emerg Med), Seidenfaden 2021 (SJTREM), Rogan 2023 (Emerg Med J), Hopman 2023 (Brain Injury), Jones 2020 (Brain Injury), Oris 2021 (J Gerontol A).

# **Table 1\. Characteristics of Included Studies**

*\[Table to be finalized with complete study characteristics including: author, year, journal, country, design, N, age, GCS range, biomarker, platform, cutoff, sampling window, CT+ prevalence, and QUADAS-2 risk of bias assessment\]*

# **Table 2\. Pooled Diagnostic Accuracy Estimates**

| Pool | k | N | Se% \[95% CI\] | Sp% \[95% CI\] | NPV% | CT avoided% | Missed /1000 |
| ----- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| GFAP+UCH-L1 | 7 | 5,903 | 96.7 \[94.7–97.9\] | 27.1 \[22.5–32.1\] | 98.75 | 24.5 | 3.1 |
|   Abbott only | 5 | 4,116 | 96.0 \[93.3–97.7\] | 27.8 \[20.8–36.2\] | 98.70 | 25.4 | 3.4 |
| S100B | 4 | 2,204 | 94.6 \[89.7–97.2\] | 25.0 \[18.2–33.4\] | 98.20 | 23.1 | 4.2 |

*Abbreviations: CT, computed tomography; NPV, negative predictive value; Se, sensitivity; Sp, specificity. All estimates from univariate random-effects models (REML) on logit scale. Abbott only \= subgroup restricted to i-STAT and Alinity platforms.*

# **Figure Legends**

**Figure 1\.** PRISMA Flow Diagram. \[To be completed\]

**Figure 2\.** Forest plots of sensitivity and specificity for GFAP+UCH-L1 (k=7, N=5,903). Pooled sensitivity was 96.7% (95% CI, 94.7%–97.9%) with I²=0%. Pooled specificity was 27.1% (95% CI, 22.5%–32.1%) with I²=93%.

**Figure 3\.** Forest plots of sensitivity and specificity for S100B (k=4, N=2,204). Pooled sensitivity was 94.6% (95% CI, 89.7%–97.2%) with I²=0%. Pooled specificity was 25.0% (95% CI, 18.2%–33.4%) with I²=93%.

**Figure 4\.** Summary receiver operating characteristic (SROC) curves comparing GFAP+UCH-L1 and S100B. Circle size is proportional to study sample size. Dashed ellipses represent 95% confidence regions.

**Figure 5\.** Clinical summary panel comparing GFAP+UCH-L1 and S100B: (A) sensitivity and specificity, (B) negative predictive value, (C) CT scans avoided and missed lesions, (D) summary table.

**eFigure 1\.** Funnel plots for assessment of publication bias for GFAP+UCH-L1 and S100B (sensitivity and specificity).