##############
# Baseline characteristics by region and arm
#
# Joshua Schwab
# jschwab77@berkeley.edu
#---------------------------

rm(list = ls())
load("outputs-withIntOnly.RData")
outputs <- outputs[resident_0 & adult_0]
nrow(outputs)

outputs[age_0 %in% 15:20, age_strata := "15-20"]
outputs[age_0 %in% 21:49, age_strata := "21-49"]
outputs[age_0 >= 50, age_strata := ">=50"]
outputs[(!is.na(htn_self_dx_0) | !is.na(htn_self_txt_0)) & !is.na(htn_uncontrol_0), htn_prev := htn_0] 
x <- expression(sex_0 == T,
                age_strata == "15-20",
                age_strata == "21-49",
                age_strata == ">=50",
                marital_0 == 1,
                marital_0 == 2,
                marital_0 == 2 & !(polygamy_0 %in% T),
                marital_0 == 3 | marital_0 == 4 | marital_0 == 5,
                is.na(marital_0),
                polygamy_0,
                is.na(polygamy_0),
                !edu_primary_0 & !edu_secondary_plus_0,
                edu_primary_0,
                edu_secondary_plus_0,
                is.na(education_0),
                formal_hi_occup_0,
                informal_hi_occup_0,
                informal_low_occup_0,
                !formal_hi_occup_0 & !informal_hi_occup_0 & !informal_low_occup_0 & occupation_0 != 18 & occupation_0 != 19,
                occupation_0 == 18 | occupation_0 == 19,
                is.na(occupation_0),
                wealth_0 == 0,
                wealth_0 == 1,
                wealth_0 == 2,
                wealth_0 == 3,
                wealth_0 == 4,
                is.na(wealth_0),
                newstable_0 == "stable",
                is.na(newstable_0),
                hiv_0,
                is.na(hiv_0),
                htn_prev,
                is.na(htn_prev)
)

GetRow <- function(x, dt) {
  if (any(grepl("htn_prev", x))) {
    dt <- dt[age_0 >= 30]
  }
  eval.x <- eval(x, dt)
  return(c(sum(eval.x, na.rm=T), mean(eval.x, na.rm = T)))
}

GetTable <- function(x, dt, col.str) {
  y <- sapply(x, GetRow, dt)
  n <- nrow(dt)
  d.out <- data.frame(c(paste0("N = ", n), paste0(format(y[1, ], big.mark=","), " (", format(round(y[2, ] * 100, 1), nsmall = 1), ")")))
  names(d.out) <- col.str
  
  return(d.out)
}

by.intervention <- cbind(GetTable(x, outputs[intervention == T], "int"), GetTable(x, outputs[intervention == F], "control"), GetTable(x, outputs, "all"))
by.region <- cbind(GetTable(x, outputs[region_name == "Western Uganda"], "W Ug"), GetTable(x, outputs[region_name == "Eastern Uganda"], "E Ug"), GetTable(x, outputs[region_name == "Kenya"], "Kenya"), GetTable(x, outputs, "all"))

rownames(by.region) <- rownames(by.intervention) <- c("", as.character(x))
filestr <- "Analysis Outputs/baseline chars.xlsx"
write.xlsx(cbind(by.intervention, by.region), file = filestr, row.names=T)

