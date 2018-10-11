##############
# Open cohort patient flow
#
# Joshua Schwab
# jschwab77@berkeley.edu
#---------------------------

rm(list = ls())
sink("Analysis Outputs/open cohort Ns.txt")
for (intervention1 in c("ALL", T, F)) {
  load("outputs-bl-with kids and move and dead.RData")
  if (intervention1 != "ALL") outputs <- outputs[intervention == intervention1]
  outputs <- outputs[!(is.na(censusdate) & tr_0==0 & chc_0==0)]
  baseline.enumerated <- nrow(outputs) 
  missing.age <- nrow(outputs[is.na(age_0)]) 
  under.15 <- nrow(outputs[age_0 < 15]) 
  outputs <- outputs[age_0 >= 15]
  baseline.adults <- nrow(outputs) 
  non.resident <- nrow(outputs[resident_0 == F | move_0 == T])
  outputs <- outputs[resident_0 == T & move_0 == F]
  died.by.start.of.chc <- nrow(outputs[dead_0 == T]) 
  outputs <- outputs[dead_0 == F] 
  baseline.adult.residents1 <- nrow(outputs)
  
  load("outputs-withIntOnly.RData") #does not include <12, move_0, dead_0
  if (intervention1 != "ALL") outputs <- outputs[intervention == intervention1]
  
  baseline.adult.residents <- nrow(outputs[age_0 >= 15 & resident_0 %in% T]) 
  stopifnot(baseline.adult.residents == baseline.adult.residents1)
  died <- nrow(outputs[age_0 >= 15 & resident_0 %in% T & dead_3 %in% T])
  outmigrated <- nrow(outputs[age_0 >= 15 & resident_0 %in% T & !(dead_3 %in% T) & (outmigrate_3 %in% T | !(resident_3 %in% T))]) 
  turned.15.years.old <- nrow(outputs[age_0 >= 12 & age_0 <= 14 & !(dead_3 %in% T) & !(outmigrate_3 %in% T) & inmigrant_3 == F & resident_3 %in% T]) 
  inmigrated <- nrow(outputs[age_0 >= 12 & !(dead_3 %in% T) & !(outmigrate_3 %in% T) & inmigrant_3 == T & resident_3 %in% T]) 
  final <- nrow(outputs[resident_3 %in% T & age_3 >= 15 & !(dead_3 %in% T) & !(outmigrate_3 %in% T)]) 
  
  results <- data.frame(baseline.enumerated, missing.age, under.15, baseline.adults, non.resident, died.by.start.of.chc, baseline.adult.residents1, died, outmigrated, turned.15.years.old, inmigrated, final)
  results <- data.frame(prettyNum(results, big.mark=","))
  names(results) <- paste0("intervention = ", intervention1)
  print(results)
}
sink()


