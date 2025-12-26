# METHODS

## DATA SOURCE AND STUDY POPULATION

This analysis utilized data from the National Health and Nutrition Examination Survey (NHANES), a nationally representative cross-sectional survey of non-institutionalized U.S. civilians conducted by the National Center for Health Statistics.^1^ NHANES employs a complex, multistage probability sampling design with oversampling of certain demographic groups to produce reliable estimates for these subpopulations.

Data from five survey cycles were included: 2005-2006, 2007-2008, 2009-2010, 2015-2016, and 2017-2020 (pre-pandemic). Cycles 2011-2014 were not incorporated due to incomplete availability of harmonized high-sensitivity C-reactive protein (hs-CRP) files for those years; their inclusion is planned as a robustness check in future analyses. The 2017-2020 cycle combined data from 2017-2018 and 2019-March 2020 (pre-pandemic portion only) per NCHS guidelines.

Participants were eligible if they were aged ≥20 years, had fasting lipid measurements (required for LDL-C calculation), had measured hs-CRP, and had prescription medication data available for statin use ascertainment. Participants with hs-CRP >10 mg/L were excluded to remove acute inflammation that may reflect transient infectious or inflammatory processes rather than chronic cardiovascular risk.^2^ Sample sizes by survey cycle are provided in Supplemental Table 1.

## VARIABLE DEFINITIONS

**Statin use.** Current statin use was identified from NHANES prescription medication files (RXQ_RX). Participants were classified as statin users if they reported taking any of the following medications in the past 30 days: atorvastatin, simvastatin, rosuvastatin, pravastatin, lovastatin, fluvastatin, or pitavastatin. Medications were identified by generic drug name matching and confirmed using the Multum Lexicon therapeutic classification code 358 (HMG-CoA reductase inhibitors).

**Low-density lipoprotein cholesterol (LDL-C).** LDL-C was calculated using the Friedewald equation: LDL-C = total cholesterol − HDL-C − (triglycerides/5), expressed in mg/dL.^3^ Calculations were restricted to participants with triglycerides <400 mg/dL, as the Friedewald equation becomes unreliable at higher triglyceride levels. Direct LDL-C measurements were used as a secondary fallback when calculated LDL-C was unavailable and triglycerides were <400 mg/dL.

**High-sensitivity C-reactive protein (hs-CRP).** hs-CRP was measured in serum using latex-enhanced nephelometry at NHANES mobile examination centers.^4^ Values are expressed in mg/L. hs-CRP ≥2 mg/L was used as the primary threshold for elevated inflammatory risk based on established cardiovascular risk stratification guidelines.^5^

**Diabetes.** Diabetes was defined as glycated hemoglobin (HbA1c) ≥6.5%, fasting plasma glucose ≥126 mg/dL, self-reported physician diagnosis of diabetes, or current use of insulin or oral hypoglycemic medications.^6^

**Prediabetes.** Prediabetes was defined as HbA1c 5.7-6.4% or fasting plasma glucose 100-125 mg/dL in the absence of diabetes.^6^

**Glycemic status.** Participants were categorized into three mutually exclusive groups: diabetes, prediabetes (excluding diabetes), and normal glycemia.

**Hypertension.** Hypertension was defined as systolic blood pressure ≥130 mm Hg, diastolic blood pressure ≥80 mm Hg, or self-reported current use of antihypertensive medication.^7^ Blood pressure was measured as the average of up to three readings obtained by trained technicians using a mercury sphygmomanometer.

**Obesity.** Obesity was defined as body mass index (BMI) ≥30 kg/m², calculated from measured weight and height.

**Current smoking.** Current smoking was defined as having smoked ≥100 cigarettes in lifetime and currently smoking cigarettes every day or some days.

## OUTCOMES

**Primary outcome.** Residual inflammatory risk (RIR) was defined as the co-occurrence of: (1) current statin use, (2) LDL-C <70 mg/dL, and (3) hs-CRP ≥2 mg/L. This definition captures patients who have achieved guideline-directed lipid control but maintain elevated systemic inflammation.

**Sensitivity definitions.** Sensitivity analyses were conducted using alternative hs-CRP thresholds (≥1.5, ≥2.5, ≥3, ≥4, ≥5 mg/L) and alternative LDL-C thresholds (<55, <60, <65, <80, <100 mg/dL) to assess the robustness of prevalence estimates.

## STRATIFICATION VARIABLES

To explore patterns of inflammatory risk across clinical subgroups, participants were stratified by:

**LDL-C and statin status (6 groups).** Participants were cross-classified into six mutually exclusive groups based on LDL-C category (≤70, 70-130, >130 mg/dL) and statin use status (user, non-user). This stratification allows comparison of inflammatory burden across the full spectrum of LDL-C levels in both statin-treated and untreated populations.

**Glycemic status.** Normal glycemia, prediabetes, and diabetes.

**Sex.** Male and female, based on self-reported sex.

## STATISTICAL ANALYSIS

All analyses incorporated NHANES complex survey design elements including primary sampling units (SDMVPSU), strata (SDMVSTRA), and sampling weights to produce nationally representative estimates.^1^ Fasting subsample weights (WTSAF2YR for cycles 2005-2016; WTSAFPRP for 2017-2020) were used given the requirement for fasting lipid measurements. Weights were divided by the number of study cycles (n=5) to appropriately combine data across multiple cycles per NCHS guidelines. The lonely primary sampling unit adjustment was applied to prevent variance estimation errors in strata with single PSUs.

**Primary analysis.** Survey-weighted prevalence of RIR was estimated among statin users with LDL-C <70 mg/dL with 95% confidence intervals. Survey-weighted logistic regression was used to identify demographic and clinical predictors of RIR. The regression model included age, sex, race/ethnicity, BMI, diabetes, hypertension, current smoking, and triglycerides as covariates. A parsimonious modeling approach was employed to avoid overfitting given the sample size.

**Exploratory analyses.** The following analyses were pre-specified as exploratory:
1. Survey-weighted prevalence of hs-CRP ≥2 mg/L across the six LDL×statin groups
2. RIR prevalence stratified by glycemic status among statin users with LDL-C <70 mg/dL
3. Sex-stratified prevalence of hs-CRP ≥2 mg/L across LDL×statin groups
4. Sex-specific logistic regression models for RIR predictors

These exploratory analyses are hypothesis-generating and were not powered for formal interaction testing.

**Characteristics comparison.** Baseline characteristics of participants with and without RIR were compared using survey-weighted means for continuous variables and survey-weighted proportions for categorical variables. P-values were calculated using survey-adjusted t-tests for continuous variables and chi-square tests for categorical variables.

**Power considerations.** The primary analytic cohort included 577 statin users with LDL-C <70 mg/dL, of whom 141 met RIR criteria. Regression analyses are limited by this sample size and are underpowered for detecting small effect sizes.

**Software.** All analyses were performed using R version 4.3 (R Foundation for Statistical Computing, Vienna, Austria). The survey package was used for survey-weighted analyses.^8^ Statistical significance was defined as two-sided P<0.05.

## ETHICAL CONSIDERATIONS

NHANES protocols were approved by the NCHS Research Ethics Review Board, and all participants provided written informed consent. This secondary analysis of publicly available, de-identified data is exempt from institutional review board review per 45 CFR 46.104.

---

## REFERENCES (Methods section)

1. Centers for Disease Control and Prevention. National Health and Nutrition Examination Survey. https://www.cdc.gov/nchs/nhanes/. Accessed December 2025.

2. Pearson TA, Mensah GA, Alexander RW, et al. Markers of inflammation and cardiovascular disease: application to clinical and public health practice: A statement for healthcare professionals from the Centers for Disease Control and Prevention and the American Heart Association. Circulation. 2003;107(3):499-511.

3. Friedewald WT, Levy RI, Fredrickson DS. Estimation of the concentration of low-density lipoprotein cholesterol in plasma, without use of the preparative ultracentrifuge. Clin Chem. 1972;18(6):499-502.

4. Rifai N, Ridker PM. High-sensitivity C-reactive protein: a novel and promising marker of coronary heart disease. Clin Chem. 2001;47(3):403-411.

5. Ridker PM. Clinical application of C-reactive protein for cardiovascular disease detection and prevention. Circulation. 2003;107(3):363-369.

6. American Diabetes Association. Classification and diagnosis of diabetes: Standards of Medical Care in Diabetes—2024. Diabetes Care. 2024;47(Suppl 1):S20-S42.

7. Whelton PK, Carey RM, Aronow WS, et al. 2017 ACC/AHA/AAPA/ABC/ACPM/AGS/APhA/ASH/ASPC/NMA/PCNA guideline for the prevention, detection, evaluation, and management of high blood pressure in adults. J Am Coll Cardiol. 2018;71(19):e127-e248.

8. Lumley T. Analysis of complex survey samples. J Stat Softw. 2004;9(8):1-19.

