##############
# Suppression and ART start
#
# Joshua Schwab
# jschwab77@berkeley.edu
#---------------------------

set.seed(42642) #for reproducibility with SuperLearner

library(ltmle)
library(data.table)
library(survival)
library(SuperLearner)

rm(list = ls())
source('Utils.R')
source('Stage2_Functions.R')
source('Adapt_Functions.R')
source('Stage1_Functions.R') 

SL.library <- list(c("SL.glm", "screen.corRank10"), c("SL.glm", "screen.corP"), c("SL.gam", "screen.corRank10"), c("SL.gam", "screen.corP"), c("SL.glm.interaction", "screen.corRank5"), c("SL.mean", "All"))

community.set <- c("Nyamrisra", "Nyatoto", "Kitare", "Ogongo", "Kisegi", "Magunga", 
                   "Sena", "Tom Mboya", "Bware", "Ongo", "Othoro", "Sibuoche", "Nankoma", 
                   "Nsiinze", "Kamuge", "Kiyunga", "Bugono", "Muyembe", "Kiyeyi", 
                   "Merikit", "Kadama", "Kameke", "Bugamba", "Nsiika", "Mitooma", 
                   "Rugazi", "Kitwe", "Rubaare", "Ruhoko", "Rwashamaire", "Kazo", 
                   "Nyamuyanja")

CreatePopulations <- function() {
  AddPop <- function(cascade = "ANY", eart_vl = "ANY", cd4 = "ANY", incorp.self) {
    pop <<- rbind(pop, data.table(cascade, eart_vl, cd4, incorp.self))  
  }
  pop <- NULL
  if (ART) {
    for (incorp.self in c(T, F)) {
      AddPop(eart_vl = F, incorp.self = incorp.self) 
      for (cd4 in 2:4) {
        AddPop(cd4 = cd4, eart_vl = F, incorp.self = incorp.self)
      }
    }
  } else {
    for (cascade in 1:4) {
      AddPop(cascade = cascade, incorp.self = F)
    }
  }
  setkeyv(pop, names(pop))
  return(pop)
}

LoadInputs <- function() {
  load("outputs-withIntOnly.RData")
  return(outputs)
}

GetDT <- function(inputs, cur.outcome, cur.pop) {
  dt <- inputs
  dt[hiv_0 & is.na(hivtest_date_0) & !is.na(VLdate_0), hivtest_date_0 := VLdate_0]
  dt[dead_3 & is.na(death_date_3), death_date_3 := HalfwayDate(hivtest_date_0, trk_end_3)] 
  dt[outmigrate_3 & is.na(outmigrate_3), outmigrate_date_3 := HalfwayDate(hivtest_date_0, trk_end_3)] 
  
  index <- dt[, resident_0 & age_0 >= 15 & hiv_0 %in% TRUE]
  if (ART) {
    index <- index & dt[, !is.na(hivtest_date_0)] 
  }
  if (cur.pop$cascade != "ANY") {
    if (cur.pop$incorp.self) {
      if (cur.pop$cascade == 1) {
        index <- index & dt[, pdx_vl_0 %in% 0 & pdx_self_0 %in% 0]
      } else if (cur.pop$cascade == 2) {
        index <- index & dt[, (pdx_vl_0 %in% 1 | pdx_self_0 %in% 1) & eart_vl_0 %in% 0]
      } else if (cur.pop$cascade == 3) {
        index <- index & dt[, (pdx_vl_0 %in% 1 | pdx_self_0 %in% 1) & (eart_vl_0 %in% 1 | eart_self_0 %in% 1) & supp_0 %in% 1]
      } else if (cur.pop$cascade == 4) {
        index <- index & dt[, (pdx_vl_0 %in% 1 | pdx_self_0 %in% 1) & (eart_vl_0 %in% 1 | eart_self_0 %in% 1) & supp_0 %in% 0]
      } else {
        stop("not expected")
      }
    } else {
      if (cur.pop$cascade == 1) {
        index <- index & dt[, pdx_vl_0 %in% 0]
      } else if (cur.pop$cascade == 2) {
        index <- index & dt[, pdx_vl_0 %in% 1 & eart_vl_0 %in% 0]
      } else if (cur.pop$cascade == 3) {
        index <- index & dt[, eart_vl_0 %in% 1 & supp_0 %in% 1]
      } else if (cur.pop$cascade == 4) {
        index <- index & dt[, eart_vl_0 %in% 1 & supp_0 %in% 0]
      } else {
        stop("not expected")
      }
    }
  }
  if (cur.pop$cd4 != "ANY") {
    if (cur.pop$cd4 == 2) {
      index <- index & dt[, !is.na(cd4_0) & cd4_0 < 350]
    } else if (cur.pop$cd4 == 3) {
      index <- index & dt[, !is.na(cd4_0) & cd4_0 >= 350 & cd4_0 < 500]
    } else if (cur.pop$cd4 == 4) {
      index <- index & dt[, !is.na(cd4_0) & cd4_0 >= 500]
    } else {
      stop("not expected")
    }
  }
  if (cur.pop$eart_vl != "ANY") {
    if (cur.pop$incorp.self) {
      index <- index & dt[, eart_vl_0 == as.logical(cur.pop$eart_vl) & eart_self_0 == as.logical(cur.pop$eart_vl)]
    } else {
      index <- index & dt[, eart_vl_0 == as.logical(cur.pop$eart_vl)]
    }
  }
  Y <- GetY(inputs, cur.outcome)
  if (ART) {
    dt$time <- Y$time
    dt$Y <- Y$Y
  } else {
    dt$Y <- Y
  }
  dt[, id := as.integer(factor(community_name))]
  dt <- dt[index]
  return(dt)
}

GetDataFrame <- function(dt) {
  if (ART) {
    df <- data.frame(time = dt$time, Y = dt$Y) 
  } else {
    df <- data.frame(get.X(as.data.frame(dt), analysis = "HIV"), C = BinaryToCensoring(is.censored = is.na(dt$Y)), Y = dt$Y)
  }
  return(df)
}

GetY <- function(inputs, cur.outcome) {
  if (cur.outcome == "supp") {
    inputs[dead_3 %in% FALSE & outmigrate_3 %in% FALSE, Y := supp_3]
    return(as.integer(inputs$Y))
  } else if (cur.outcome %in% c("ART6", "ART12", "ART24", "ART36")) {
    event.time <- inputs[, art_date]
    censor.time <- inputs[, pmin(death_date_3, outmigrate_date_3, trk_end_3, na.rm = T)]
    event <- DateBefore(event.time, censor.time, strict = F)
    time <- as.numeric(pmin(event.time, censor.time, na.rm = T) - inputs[, hivtest_date_0], units = "days")
    return(list(time=time, Y=event))
  } else {
    stop("not expected")
  }
}

Est <- function(cur.outcome, cur.pop) {
  dt <- GetDT(inputs, cur.outcome, cur.pop)
  estimates <- data.table(community_name = community.set, Y = NA_real_)
  for (comm in community.set) {
    df <- GetDataFrame(dt[community_name == comm])
    stopifnot(nrow(df) > 1 & sum(!is.na(df$Y)) > 1)
    if (ART) {
      stopifnot(!anyNA(df$Y))
      stopifnot(!anyNA(df$time))
      months <- as.integer(substr(cur.outcome, start = 4, stop = nchar(cur.outcome)))
      days <- 365.25 / 12 * months
      
      fit <- survfit(Surv(time, Y) ~ 1, data = df)
      est <- 1 - summary(fit, min(days, max(df$time)))$surv
    } else {
      est <- ltmle(df, Anodes = NULL, abar = NULL, Cnodes = "C", Ynodes = "Y", estimate.time = F, variance.method = "ic", SL.library = SL.library)$estimates["tmle"]
    }
    estimates[community_name == comm, Y := est]
  }
  return(estimates)
}

Get1 <- function(z) {
  unique.z <- unique(z)
  stopifnot(length(unique.z) == 1)
  return(unique.z)
}

inputs <- LoadInputs()
for (ART in c(T, F)) {
  if (ART) {
    outcomes <- paste0("ART", c(6, 12, 24, 36))
  } else {
    outcomes <- "supp"
  }
  
  populations <- CreatePopulations()
  results.stage2 <- data.table()
  for (outcome.index in 1:length(outcomes)) {
    cur.outcome1 <- outcomes[outcome.index]
    
    for (i in 1:nrow(populations)) {
      cur.pop1 <- populations[i]
      
      cat("\n\n------ outcome.index = ", outcome.index, " i = ", i, "\n")
      print(cur.pop1)
      
      dt1 <- GetDT(inputs, cur.outcome1, cur.pop1)
      estimates <- Est(cur.outcome1, cur.pop1)
      dt.summary <- dt1[, .(below25 = mean(age_0 < 25), A = Get1(as.integer(intervention)), pair = Get1(as.integer(as.character(pair))), supp = mean(supp_0, na.rm=T), nIndv = .N, U = 1, id = Get1(id)), by = "community_name"]
      
      dt.summary.merged <- merge(estimates, dt.summary, by = "community_name")
      stage2 <- Stage2(data.input=as.data.frame(dt.summary.merged), weighting = "individual", outcome = "Y", clust.adj = c("U", "below25", "supp"), do.data.adapt = T, break.match = F)
      results.stage2 <- rbind(results.stage2, cbind(cur.pop1, outcome = cur.outcome1, N = nrow(dt1), stage2))
    }
  }
  filestr <- paste("Analysis Outputs/closed cohort outputs - ART =", ART)
  write.csv(format(results.stage2, digits = 3), file=paste(filestr, "- stage2.csv"))
}
