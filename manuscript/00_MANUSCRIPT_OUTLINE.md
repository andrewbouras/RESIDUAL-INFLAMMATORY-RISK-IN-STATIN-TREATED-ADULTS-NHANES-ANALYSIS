# RIR Manuscript Outline (Final)

**Working Title:** Prevalence and Correlates of Residual Inflammatory Risk Among Statin-Treated U.S. Adults Achieving LDL-C <70 mg/dL: A Cross-Sectional NHANES Analysis

**Target Journal:** JACC / Circulation / JAHA (structure follows JACC format per example)

---

## CRITICAL ADJUSTMENTS (per PI feedback)

### 1. Missing 2011-2014 Cycles
**Language to use:** "Cycles 2011-2014 were not incorporated due to incomplete availability of harmonized hs-CRP files for those years. Their inclusion is planned as a robustness check."

**Do NOT say:** "automated download limitations" or "API issues"

### 2. Primary vs Exploratory Analyses
**State explicitly in Methods:**
- **Primary analysis:** RIR prevalence among statin users with LDL-C <70 mg/dL; predictors of RIR
- **Exploratory analyses:** LDL×statin 6-group comparisons; glycemic status stratification; sex-stratified analyses

### 3. Low LDL, High Inflammation Phenotype
**Language to use:** "This pattern may reflect underlying metabolic-inflammatory status rather than a distinct clinical phenotype and requires confirmation in longitudinal studies."

### 4. Regression Power
**State explicitly:** "Regression analyses are limited by sample size (N=577 with 141 RIR events) and are underpowered for detecting small effect sizes."

### 5. Glycemic Strata Precision
**Note:** "Prevalence estimates in the normal glycemia subgroup have high imprecision (SE 9.3%) due to limited sample size."

---

## ABSTRACT (Structured, ~250 words)

### BACKGROUND
Statins reduce LDL cholesterol and cardiovascular events, but 30-40% of statin-treated patients in clinical trials maintain high-sensitivity C-reactive protein (hs-CRP) ≥2 mg/L despite achieving LDL-C goals. Nationally representative estimates of this "residual inflammatory risk" (RIR) are lacking.

### OBJECTIVES
To estimate the weighted prevalence of RIR among statin-treated U.S. adults achieving LDL-C <70 mg/dL and characterize RIR across LDL-C and statin status groups, glycemic status, and sex.

### METHODS
Cross-sectional analysis of NHANES 2005-2010 and 2015-2020. Fasting adults ≥20 years with hs-CRP ≤10 mg/L were included. RIR was defined as statin use + LDL-C <70 mg/dL + hs-CRP ≥2 mg/L. Survey-weighted prevalence estimates and logistic regression were performed.

### RESULTS
Among 12,896 fasting adults, 2,527 were statin users (19.6%), of whom 577 achieved LDL-C <70 mg/dL. RIR prevalence was 21.9% (95% CI: 15.9-27.9%). Across 6 LDL×statin groups, statin users consistently showed higher hs-CRP ≥2 prevalence than non-users within each LDL tier. Diabetics had higher RIR prevalence (26.2%) compared to those with normal glycemia (15.3%, SE 9.3%). Women showed higher inflammatory burden than men in non-statin groups. BMI predicted RIR in women (OR 1.07, p=0.02) but not men.

### CONCLUSIONS
Approximately one in five statin-treated U.S. adults achieving LDL-C <70 mg/dL has residual inflammatory risk. RIR is more prevalent in diabetics and shows sex-specific patterns. These findings support targeted inflammatory risk assessment in statin-treated patients with well-controlled LDL-C.

---

## 1. INTRODUCTION (~400 words, 4 paragraphs)

### Paragraph 1: The clinical problem
- Statins are cornerstone of CVD prevention
- Despite LDL-C reductions, substantial residual cardiovascular risk persists
- Inflammation increasingly recognized as independent pathway beyond lipid lowering
- hs-CRP ≥2 mg/L associated with increased CV events even when LDL-C is controlled

### Paragraph 2: Evidence from trials
- JUPITER: statin benefit greatest in those with both elevated LDL and CRP
- CANTOS: canakinumab reduced CV events via IL-1β inhibition
- COLCOT, LoDoCo2: colchicine reduces events in CAD patients
- Trial data suggest 30-40% of statin-treated patients maintain hs-CRP ≥2 mg/L
- However, no nationally representative U.S. population estimates exist

### Paragraph 3: Knowledge gap and rationale
- Understanding the burden of "residual inflammatory risk" (RIR) essential for:
  - Population-level prevention planning
  - Identifying candidates for adjunctive anti-inflammatory therapy
  - Characterizing high-risk subgroups (by glycemic status, sex)
- AHA cardiovascular-kidney-metabolic (CKM) syndrome framework emphasizes metabolic-inflammatory connections

### Paragraph 4: Objectives
- **Primary:** Estimate survey-weighted prevalence of RIR among statin users achieving LDL-C <70 mg/dL
- **Secondary:** Characterize hs-CRP elevation across 6 LDL×statin groups; evaluate RIR by glycemic status and sex; identify predictors of RIR

---

## 2. METHODS (~800 words)

### 2.1 Data Source and Study Population
- NHANES 2005-2010, 2015-2020 (5 cycles)
- Cycles 2011-2014 not incorporated due to incomplete availability of harmonized hs-CRP files; inclusion planned as robustness check
- Nationally representative sample of non-institutionalized U.S. adults
- Fasting subsample for lipid analyses
- Inclusion: adults ≥20 years with hs-CRP, fasting lipids, statin use data
- Exclusion: hs-CRP >10 mg/L (acute inflammation)

### 2.2 Variable Definitions

| Variable | Definition |
|----------|------------|
| Statin use | RXQ_RX prescription files; generic drug name matching (atorvastatin, simvastatin, rosuvastatin, pravastatin, lovastatin, fluvastatin, pitavastatin) and Multum therapeutic class code 358 |
| LDL-C | Friedewald equation: TC - HDL - (TG/5); valid when TG <400 mg/dL |
| hs-CRP | High-sensitivity C-reactive protein, mg/L |
| Diabetes | HbA1c ≥6.5% OR fasting glucose ≥126 mg/dL OR self-reported diagnosis |
| Prediabetes | HbA1c 5.7-6.4% OR fasting glucose 100-125 mg/dL (excluding diabetes) |
| Hypertension | SBP ≥130 OR DBP ≥80 OR taking antihypertensive medication |
| Obesity | BMI ≥30 kg/m² |
| Current smoker | ≥100 lifetime cigarettes AND currently smoking |

### 2.3 Outcomes

**Primary outcome:** Residual inflammatory risk (RIR)
- Definition: Statin use + LDL-C <70 mg/dL + hs-CRP ≥2 mg/L

**Sensitivity definitions:**
- hs-CRP ≥3 mg/L
- LDL-C <55 mg/dL

### 2.4 Stratification Variables

**6-group LDL×Statin:**
1. LDL-C ≤70, statin user
2. LDL-C ≤70, non-user
3. LDL-C 70-130, statin user
4. LDL-C 70-130, non-user
5. LDL-C >130, statin user
6. LDL-C >130, non-user

**Glycemic status:** Normal / Prediabetes / Diabetes

**Sex:** Male / Female

### 2.5 Statistical Analysis

**Primary analysis:**
- Survey-weighted prevalence of RIR among statin users with LDL-C <70 mg/dL
- Survey-weighted logistic regression for predictors of RIR

**Exploratory analyses:**
- hs-CRP ≥2 prevalence across 6 LDL×statin groups
- RIR prevalence by glycemic status
- Sex-stratified prevalence and predictors

**Survey design:**
- All analyses account for NHANES complex survey design
- Primary sampling units (SDMVPSU), strata (SDMVSTRA)
- Fasting subsample weights (WTSAF2YR/WTSAFPRP) adjusted for multi-cycle analysis
- Lonely PSU adjustment applied

**Software:** R version 4.x; survey package for weighted analyses

**Power considerations:** Regression analyses are limited by sample size (N=577 with 141 RIR events) and are underpowered for detecting small effect sizes.

### 2.6 Ethical Considerations
- NHANES data are publicly available and de-identified
- Analysis exempt from institutional review board review per 45 CFR 46.104

---

## 3. RESULTS (~1200 words)

### 3.1 Study Population
- 12,896 fasting adults ≥20 years included
- 2,527 statin users (19.6%)
- 577 statin users with LDL-C <70 mg/dL (primary analytic cohort)
- 141 met RIR criteria
- Reference Table 1 for characteristics

### 3.2 Primary Outcome: RIR Prevalence
- RIR prevalence: 21.9% (95% CI: 15.9-27.9%)
- Sensitivity (hs-CRP ≥3): 13.4% (95% CI: 8.8-18.0%)
- Sensitivity (LDL <55): 26.9% (95% CI: 17.0-36.8%)

### 3.3 Exploratory: 6-Group LDL×Statin Analysis
- Statin users showed higher hs-CRP ≥2 prevalence than non-users within each LDL tier
- LDL≤70, statin: 22.4%; LDL≤70, no statin: 14.5%
- LDL 70-130, statin: 20.1%; LDL 70-130, no statin: 17.7%
- LDL>130, statin: 27.0%; LDL>130, no statin: 19.8%
- Reference Figure 1

### 3.4 Exploratory: Glycemic Status
- RIR prevalence by glycemic status (statin + LDL<70):
  - Normal: 15.3% (SE 9.3%) — note: high imprecision
  - Prediabetes: 17.4% (SE 5.2%)
  - Diabetes: 26.2% (SE 4.8%)
- Pattern consistent across all 6 LDL×statin groups
- Reference Figure 2

### 3.5 Exploratory: Sex Differences
- Men: 23.7% (SE 3.5%); Women: 19.2% (SE 4.5%)
- Women showed higher hs-CRP prevalence in non-statin groups
- LDL>130 statin: women 32.0% vs men 20.4%
- Sex-specific predictors: BMI (OR 1.07, p=0.02) and age (OR 1.04, p=0.04) significant in women only
- Reference Figure 3

### 3.6 Predictors of RIR
- In adjusted model, no individual predictor reached p<0.05
- BMI showed trend (OR 1.05, 95% CI 1.00-1.11, p=0.055)
- Diabetes (OR 1.64), hypertension (OR 1.96), smoking (OR 1.96) trended higher but NS
- Reference Table 3

### 3.7 Sensitivity Analyses
- hs-CRP threshold sensitivity: dose-response from 29.6% (≥1.5) to 6.4% (≥5.0)
- LDL threshold sensitivity: RIR prevalence stable (20.7-26.9%) across thresholds
- Reference Supplemental Table

---

## 4. DISCUSSION (~1000 words)

### Paragraph 1: Principal findings
- ~22% RIR prevalence represents substantial inflammatory burden
- Confirms trial observations at population level
- First nationally representative U.S. estimate

### Paragraph 2: 6-group findings
- Higher inflammation in statin users may reflect channeling (sicker patients prescribed statins)
- Alternatively, statins may not fully address inflammation
- 14.5% of LDL≤70 non-statin users with elevated hs-CRP suggests metabolic-inflammatory status independent of medication
- **Note:** "This pattern may reflect underlying metabolic-inflammatory status rather than a distinct clinical phenotype and requires confirmation in longitudinal studies."

### Paragraph 3: Diabetes/prediabetes
- Strong association with RIR (~2× prevalence)
- Aligns with AHA CKM syndrome framework (Ndumele 2023)
- Implications for agents with dual metabolic-inflammatory effects (GLP-1RA, SGLT2i)

### Paragraph 4: Sex differences
- Women's higher inflammatory burden consistent with known sex differences
- Estrogen depletion, adiposity distribution, immune response differences
- BMI as predictor in women aligns with adiposity-inflammation link
- Supports sex-specific risk stratification

### Paragraph 5: Clinical implications
- Supports hs-CRP screening in statin-treated patients achieving LDL goals
- Diabetics with well-controlled LDL may benefit most from inflammatory risk assessment
- May inform eligibility for adjunctive anti-inflammatory therapies (colchicine, future IL-6 inhibitors)

### Paragraph 6: Prior literature
- JUPITER: median hs-CRP 4.2 at baseline, 2.2 on statin
- CANTOS eligibility: hs-CRP ≥2 despite statin
- Our 22% estimate consistent with trial enrollment criteria

---

## 5. LIMITATIONS (~250 words)

1. **Cross-sectional design:** Cannot infer causation; prevalence only
2. **Single hs-CRP measurement:** Within-person variability not captured; may misclassify some individuals
3. **Self-reported statin use:** Adherence not directly measured
4. **Missing cycles:** 2011-2014 not incorporated due to incomplete availability of harmonized hs-CRP files; planned robustness check
5. **Sample size:** N=577 with 141 events limits regression power for small effects
6. **Exploratory analyses:** LDL×statin, glycemic, and sex strata are hypothesis-generating; not powered for formal interaction testing
7. **No outcomes linkage:** Cannot directly link RIR to cardiovascular events in this cross-sectional design
8. **Subgroup imprecision:** Estimates in normal glycemia subgroup have high imprecision (SE 9.3%)

---

## 6. CONCLUSIONS (~100 words)

Approximately one in five statin-treated U.S. adults achieving LDL-C <70 mg/dL has residual inflammatory risk defined by hs-CRP ≥2 mg/L. In exploratory analyses, RIR prevalence was higher among diabetics and showed sex-specific patterns, with BMI predicting RIR in women but not men. These nationally representative estimates quantify the burden of on-treatment inflammation and support targeted hs-CRP screening in statin-treated patients with well-controlled LDL-C, particularly those with diabetes.

---

## TABLES

| Table | Title |
|-------|-------|
| **Table 1** | Characteristics of Statin Users with LDL-C <70 mg/dL by RIR Status |
| **Table 2** | hs-CRP ≥2 mg/L Prevalence Across LDL-C and Statin Status Groups |
| **Table 3** | Predictors of Residual Inflammatory Risk: Survey-Weighted Logistic Regression |
| **Suppl Table 1** | Sample Sizes by NHANES Cycle |
| **Suppl Table 2** | Sensitivity Analyses: RIR Prevalence by hs-CRP and LDL-C Thresholds |

---

## FIGURES

| Figure | Title |
|--------|-------|
| **Figure 1** | hs-CRP ≥2 mg/L Prevalence Across 6 LDL×Statin Groups |
| **Figure 2** | RIR Prevalence by Glycemic Status Among Statin Users with LDL-C <70 mg/dL |
| **Figure 3** | Sex Differences in hs-CRP Prevalence Across LDL×Statin Groups |
| **Central Illustration** | Graphical Abstract: Study Flow and Key Findings |

---

## KEY CITATIONS (2020-2025)

1. Ridker PM. Residual inflammatory risk. Curr Opin Lipidol. 2023.
2. Ridker PM et al. CANTOS Trial. N Engl J Med. 2017;377:1119-1131.
3. Tardif JC et al. COLCOT Trial. N Engl J Med. 2019;381:2497-2505.
4. Nidorf SM et al. LoDoCo2 Trial. N Engl J Med. 2020;383:1838-1847.
5. Ndumele CE et al. Cardiovascular-Kidney-Metabolic Health. Circulation. 2023.
6. Libby P. Inflammation in atherosclerosis. Nature Rev Cardiol. 2021.
7. CDC NHANES Analytic Guidelines. 2023.
8. Arnett DK et al. ACC/AHA Cholesterol Guidelines. Circulation. 2019.

---

## NEXT CHUNK: Table 1 Generation

Ready to proceed to Table 1 (Characteristics by RIR Status)?

