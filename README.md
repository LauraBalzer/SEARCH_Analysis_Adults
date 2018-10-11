
This repository contains the functions to implement the Two-Stage analysis for estimating 
the effect of the SEARCH intervention on adult outcomes. 
The SEARCH statistical analysis plan is available at https://arxiv.org/abs/1808.03231.

R Code by Laura B. Balzer, Yea-Hung Chen, and Joshua Schwab 


Stage 1: Estimate the community-level outcome. 
- Cumulative HIV incidence: HIV_Master.R, HIV_Functions.R
- NCD analyses: NCD_Master.R, NCD_Functions.R
- Population-level viral suppression & UNAIDS 90/90/90: Cascade_Master.R, Cascade_Functions.R

- Suppression and ART start among closed cohort: SuppressionInClosedCohort.R
- Annual HIV incidence in the intervention arm: Estimate Incidence.R
- Cumulative HIV testing: Cumulative Testing.R

- Mortality: mortality analysis.r
- Tuberculosis analyses: tuberculosis analysis.r


Stage 2: Estimate the intervention effect, given community-level covariates, exposure, and outcomes. 
- Stage2_Functions.R
- Adapt_Functions.R


Also includes code for descriptive analyses described in the statistical analysis plan:
- baseline chars.R
- baseline prev by age and sex and region.R  
- patient flow primary endpoint.Rmd
- Open Cohort Patient Flow.R
- Annual HIV Incidence Cohorts in Intervention Arm.Rmd
- changes in Circumcision.Rmd

#——--
Helper functions used in analyses for HIV incidence, NCD, suppression/coverage: 
- Stage1_Functions 

Helper functions used in mortality/tuberculosis analyses: 
- impute dates.r
- define survival variables.r
- community-level estimates.r
- Preprocess_Functions.R
- s2.r

Helper functions used in analyses for cumulative testing, suppression/ART in the closed cohort, 
estimating annual incidence in the intervention arm
- Utils.R
