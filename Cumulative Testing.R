##############
# Cumulative HIV testing
#
# Joshua Schwab
# jschwab77@berkeley.edu
#---------------------------

rm(list=ls())
source('Utils.R')

load("outputs-withIntOnly.RData")

outputs[is.na(outmigrate_date_3) & outmigrate_3, outmigrate_date_3 := HalfwayDate(chc_start_0, chc_start_3)]
outputs[is.na(death_date_3) & dead_3, death_date_3 := HalfwayDate(chc_start_0, chc_start_3)]

results <- NULL
for (int in c(T, F)) {
  arm <- if (int) "I" else "C"
  for (t in 0:3) {
    dt <- outputs[, grep("age|chc_start_|all_hiv_|death_date_3|inmigrant_date_3|outmigrate_date_3|intervention|self_hivtest_0|resident_", names(outputs)), with = F]
    dt$age_t <- dt[[paste0("age_", t)]]
    dt$chc_start_t <- dt[[paste0("chc_start_", t)]]
    if (t %in% c(0, 3)) {
      dt$resident_t <- dt[[paste0("resident_", t)]]
    } else {
      dt$resident_t <- dt[[paste0("resident_recensus_", t)]]
    }
    dt <- dt[intervention == int & resident_t & age_t >= 15 & !DateBefore(death_date_3, chc_start_t) & !DateBefore(outmigrate_date_3, chc_start_t)]
    
    dt[, ever_test_t := F]
    for (j in 0:t) {
      dt$ever_test_t <- dt$ever_test_t | !is.na(dt[[paste0("all_hiv_", j)]])
    }
    if (t == 0) {
      results <- rbind(results, data.table(arm, t = "self", p = dt[, mean(self_hivtest_0, na.rm = T)]))
      cat("t = 0, int = ", int, " dt[!is.na(self_hivtest_0), .N] = ", dt[!is.na(self_hivtest_0), .N], "\n")
    }
    results <- rbind(results, data.table(arm, t, p = dt[, mean(ever_test_t)]))
    cat("t = ", t, " int = ", int, " dt[, .N] = ", dt[, .N], "\n")
  }
}
write.csv(results, "Analysis Outputs/cumulative testing.csv")