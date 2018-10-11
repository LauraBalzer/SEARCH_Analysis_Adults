##############
# HIV Incidence over time
#
# Joshua Schwab
# jschwab77@berkeley.edu
#---------------------------

library(sandwich)
rm(list = ls())
source('Utils.R')

CreatePopulations <- function() {
  AddPop <- function(stable, region = "any", sex = "any") {
    pop <<- rbind(pop, data.table(stable, region, sex))  
  }
  
  pop <- NULL
  for (stable in c("any", T)) {
    AddPop(stable = stable)
    for (region in c("Western Uganda", "Eastern Uganda", "Kenya", "any")) {
      AddPop(stable = stable, region = region)
      if (region == "any") {
        for (sex in c(T, F)) {
          AddPop(stable = stable, region = region, sex = sex)
        }
      }
    }
  }
  return(unique(pop))
}

LoadInputs <- function() {
  load("outputs-withIntDeath.RData")
  inputs <- outputs[!is.na(sex_0) & intervention]
  
  #if all_hiv_t is TRUE but hivtest_date_t is missing, set to earliest evidence of hiv
  inputs[, trk_end_0 := ext_trk_end_0]
  for (tt in 0:3) {
    inputs[get(paste0("hivtest_date_", tt)) < get(paste0("chc_start_", tt)),  paste0("hivtest_date_", tt) := get(paste0("chc_start_", tt))]
    inputs[get(paste0("hivtest_date_", tt)) > get(paste0("trk_end_", tt)),  paste0("hivtest_date_", tt) := get(paste0("trk_end_", tt))]
    
    inputs[get(paste0("hivtest_", tt)) == "Positive", hivtest_date_pos := get(paste0("hivtest_date_", tt))]
    
    inputs[get(paste0("hiv_preCHC_", tt)) %in% T, paste0("hivtest_date_", tt) := pmin(hivtest_date_pos, min.prior.bc.date, tb_hiv_date_moh, get(paste0("chc_start_", tt)), na.rm=T)]
    inputs[get(paste0("all_hiv_", tt)) %in% T & is.na(get(paste0("hivtest_date_", tt))), paste0("hivtest_date_", tt) := get(paste0("chc_start_", tt))]
  }
  
  inputs[, dead_1 := dead_s_1]
  inputs[, dead_2 := dead_1 | dead_s_2]
  inputs[, dead_3 := dead_2 | dead_3]
  
  inputs[searchid %in% c("27109305001-1", "36317578051-10", "42135601067-2"), hivtest_date_1 := pmin(hivtest_date_1, hivtest_date_2)] #correct errors in hivtest_date_1 for 3 searchids
  return(inputs)
}

GetDT <- function(inputs, tt, cur.pop) {
  stopifnot(tt %in% 0:2)
  next.t <- tt + 1
  dt <- inputs
  Y.list <- GetY(inputs, tt)
  Y.list_next <- GetY(inputs, next.t)
  dt$Y_t <- Y.list$Y
  dt$Y_next <- Y.list_next$Y
  dt$hivtest_date_t <- Y.list$hivtest_date
  dt$hivtest_date_next <- Y.list_next$hivtest_date
  
  dt[, num.years := as.numeric(hivtest_date_next - hivtest_date_t, units = "days") / 365.25]
  dt[Y_next == 1, pt := num.years / 2]
  dt[Y_next == 0, pt := num.years]
  
  age <- dt[[paste0("age_", tt)]]
  resident <- dt[[paste0("resident_", tt)]]
  index <- age >= 15 & resident
  if (cur.pop$stable == "TRUE") {
    index <- index & dt[, stable_0 %in% T]
  }
  if (cur.pop$region != "any") {
    index <- index & dt[, region_name == cur.pop$region]
  }
  if (cur.pop$sex != "any") {
    index <- index & dt[, sex_0 == cur.pop$sex]
  }
  dt <- dt[index]
  dt[, current.t := tt]
  return(dt)
}

GetAdjustmentSet <- function(dt) {
  tt <- unique(dt$current.t)
  stopifnot(length(tt) == 1 && tt %in% 0:2)
  adj.set <- dt[, c("mobile_0", "sex_0", paste0("age_", tt)), with = F]
  setnames(adj.set, paste0("age_", tt), "age_t")
  return(adj.set)
}

GetDataFrame <- function(dt, A) {
  index <- dt[, Y_t %in% 0 & !is.na(Y_next)]
  df <- data.frame(GetAdjustmentSet(dt), A, Y_next = dt[, Y_next])[index, ]
  return(list(df = df, searchid = dt$searchid[index], pt = dt$pt[index]))
}

GetY <- function(inputs, tt) {
  Y <- as.integer(inputs[[paste0("all_hiv_", tt)]])
  hivtest_date <- inputs[[paste0("hivtest_date_", tt)]]
  stopifnot(!anyNA(hivtest_date[Y %in% F]))
  Y[inputs[[paste0("dead_", tt)]] %in% T] <- NA_integer_
  hivtest_date[is.na(Y)] <- NA_Date
  return(list(Y = Y, hivtest_date = hivtest_date))
}

EstIncidenceRate <- function(cur.pop) {
  inc.rate <- numeric(3)
  for (t in 0:2) {
    dt <- GetDT(inputs, t, cur.pop)
    inc.rate[t + 1] <- dt[intervention & !is.na(Y_next) & Y_t == 0, sum(Y_next) / sum(pt)]
  }
  return(inc.rate)
}

EstIncidenceChange <- function(cur.pop) {
  df0.list <- GetDataFrame(GetDT(inputs, tt = 0, cur.pop), A = 0)
  dfA.list <- GetDataFrame(GetDT(inputs, tt = 2, cur.pop), A = 1)
  
  df.pooled <- rbind(df0.list$df, dfA.list$df)
  searchid <- c(as.character(df0.list$searchid), as.character(dfA.list$searchid))
  offset1 <- c(log(df0.list$pt), log(dfA.list$pt))
  ind.poisson.rate.results <- EstIncidenceChangePoisson(df.pooled, id = searchid, offset1 = offset1)
  
  return(ind.poisson.rate.results)
}

EstIncidenceChangePoisson <- function(df, id, offset1, adjust) {
  stopifnot(setequal(names(df), c("mobile_0", "sex_0", "age_t", "A", "Y_next")))
  m <- glm(Y_next ~ mobile_0 + sex_0 + age_t + A + offset(offset1), family=poisson(link="log"), data=data.frame(df, offset1))
  
  est <- coef(m)["A"] 
  v <- vcovCL(m, cluster=id) #same as gee with "exchangeable"
  s <- sqrt(diag(v))["A"]
  
  stopifnot(uniqueN(id) > 100) #use t dist if smaller
  x <-  qnorm(0.975) * s
  CI <- c("2.5%" = est - x, "97.5%" = est + x)
  pvalue <- 2 * pnorm(-abs(est / s))
  return(c(exp(est), exp(CI), pvalue))
}

GetDescriptiveStats <- function(cur.pop) {
  results <- NULL
  for (t in 0:2) {
    dt <- GetDT(inputs, t, cur.pop)
    uninfected <- dt[Y_t == 0, .N]
    df.list <- GetDataFrame(dt, A=NA) #only includes Y_t == 0 and !is.na(Y_t+1)
    pt <- sum(df.list$pt) 
    outcome.measured <- nrow(df.list$df)
    results <- rbind(results, data.table(cur.pop, t, uninfected, pt, outcome.measured))
  }
  return(results)
}

inputs <- LoadInputs()
populations <- CreatePopulations()

num.populations <- nrow(populations)
results.change <- array(dim=c(num.populations, 4))
results.rate <- array(dim=c(num.populations, 3))
colnames(results.change) <- c("RR", "CIl", "CIh", "p")
colnames(results.rate) <- paste0("t = ", 0:2)
results.stats <- NULL
for (i in 1:num.populations) {
  cat("i = ", i, "\n")
  cur.pop <- populations[i]
  results.rate[i, ] <- EstIncidenceRate(cur.pop)
  results.change[i, ] <- EstIncidenceChange(cur.pop)
  results.stats <- rbind(results.stats, GetDescriptiveStats(cur.pop))
}

path <- "Analysis Outputs/"
write.csv(format(data.frame(populations, results.change), digits = 3), file=paste0(path, "incidence results-change.csv"), row.names = F)
write.csv(format(data.frame(populations, results.rate, check.names = F), digits = 3), file=paste0(path, "incidence results-rate.csv"), row.names = F)
write.csv(results.stats, file=paste0(path, "incidence results-stats.csv"), row.names = F)
