##############
# Baseline prevalence
#
# Joshua Schwab
# jschwab77@berkeley.edu
#---------------------------

library(dplyr)
rm(list = ls())
load("outputs-withIntOnly.RData")
outputs[age_0 < 15, age_strata := "12-14"]
outputs[age_0 %in% 15:24, age_strata := "15-24"]
outputs[age_0 %in% 25:34, age_strata := "25-34"]
outputs[age_0 %in% 35:44, age_strata := "35-44"]
outputs[age_0 %in% 45:54, age_strata := "45-54"]
outputs[age_0 >= 55, age_strata := "55+"]
outputs[!is.na(sex_0), sex.str := if_else(sex_0, "male", "female")]
out <- outputs[!is.na(sex_0) & age_0 >= 15 & !is.na(hiv_0) & resident_0, .(baseline_prev = mean(hiv_0), num = sum(hiv_0), denom = .N), keyby=c("region_name", "age_strata", "sex.str")]
write.csv(out, file = "Analysis Outputs/baseline prev by age and sex and region.csv", row.names = F)
