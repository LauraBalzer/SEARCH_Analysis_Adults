#==========================================
# HIV FUNCTIONS
# R code to complete the analyses for the HIV primary & secondary outcomes described
#	in detail in the SEARCH Analysis Plan

# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Lead Statistician for SEARCH
#==========================================



Run.HIV<- function(subgroup, SL.library){
  
  
  data <- preprocess.HIV()
  print( dim(data))
  
  # Stage 1 - estimate HIV incidence within each community & output aggregated data
  data  <- Stage1.HIV(data,  subgroup=subgroup,
                      SL.library=SL.library)
  # save covariates
  W<- data$W
  data <- data$data.clust
  dim(data)
  
  
  # Stage 2- estimation of the intervention effect on HIV incidence
  
  # candidate adjustment variables 
  if(subgroup=='Kenya' | subgroup=='EU' | subgroup=='SWU'){
    clust.adj <- 'U'
  } else{
    clust.adj <- c('U', 'hiv_prev_0', 'mc_cover_0')
  } 
  
  # primary is arithmetic risk ratio  
  goal <- 'aRR'
  # primary weights  clusters (communities)  equally
  weighting <-  'clust'

  # primary stage1 outcome Yc is unadjusted incidence
  # primary stage2 uses adaptive prespecification

  data.temp <- reformat.stage2(data, outcome='inc.U', indv='oic.size')
  # primary analysis
  unadj.inc_adj.effect <- Stage2(goal=goal, weighting=weighting, data.input= data.temp, 
                                 outcome='Y', clust.adj=clust.adj,  do.data.adapt=T,  
                                 break.match= F)
  # secondary - no adjustment in Stage2
  unadj.inc_unadj.effect <- Stage2(goal=goal, weighting=weighting, data.input= data.temp, 
                                   outcome='Y', clust.adj=NULL, do.data.adapt=F,  break.match= F)
  
  
  # secondary stage1 outcome Yc is adjusted incidence
  data.temp <- reformat.stage2(data, outcome='inc.A', indv='oic.size')
  
  adj.inc_adj.effect <- Stage2(goal=goal, weighting=weighting, data.input= data.temp, 
                               outcome='Y', clust.adj=clust.adj,  do.data.adapt=T,  break.match= F)
  adj.inc_unadj.effect <- Stage2(goal=goal, weighting=weighting, data.input= data.temp, 
                                 outcome='Y', clust.adj=NULL,  do.data.adapt=F,  break.match= F)
  
  # crude incidence rate
  data.temp <- reformat.stage2(data, outcome='IR', indv='PY.tot')
  weighting <-  'indv'
  
  crude.inc.rate <- Stage2(goal=goal, weighting=weighting, data.input= data.temp, 
                           outcome='Y', clust.adj=NULL, do.data.adapt=F,  break.match= F)
  
  
  est <- data.frame(rbind(unadj.inc_adj.effect, unadj.inc_unadj.effect, 
                          adj.inc_adj.effect, adj.inc_unadj.effect,
                          crude.inc.rate))
  
  
  rownames(est) <-  c('Unadj-Incidence; Adj-Effect',
                      'Unadj-Incidence; Unadj-Effect',
                      'Adj-Incidence; Adj-effect',
                      'Adj-Incidence; Unadj-Effect',
                      'Crude Incidence Rate')
  
  
  rownames(data) <- data$community_name
  data<- subset(data, select=-c(U,community_name, id))
  data[data$region_name=='Western Uganda','region_name'] <- 3
  data[data$region_name=='Eastern Uganda', 'region_name'] <- 2
  data[data$region_name=='Kenya', 'region_name'] <- 1
  colnames(data)[1]<- 'region'
  
  data$region <- as.numeric(data$region)
  RETURN <- list(W=W, data=data, est=est)
  
  RETURN
}



# preprocess: 
preprocess.HIV <- function(){
  
  load("outputs-withIntOnly.RData")
  print( dim(outputs) )
  data.input <- outputs
  
  # Exclude anyone that was flagged as an SEARCH-id related error
  data.input <- subset(data.input, !(data_flag | dead_0 | move_0) )
  
  # outputs-withIntOnly.Rdata does NOT have community_number
  # adding now
  id <- rep(NA, nrow(data.input))
  comm <- unique(data.input$community_name)
  for(j in 1:length(comm)){
    
    these.units <- data.input$community_name==comm[j]
    id[these.units] <- j
    
  }	
  
  data.input <- cbind(data.input, 
                      A= as.numeric(as.logical(data.input$intervention)),
                      id=id
                      )
  
  # transform pairs to be numeric
  data.input$pair <- as.numeric(as.character(data.input$pair))
  
  # if haven't added a dummy variable for unadjusted
  if( sum(grep('U', colnames(data.input)))==0){
    data.input <- cbind(U=1, data.input) 
  }
  
  print('***preprocessing done***')
  
  data.input
}



# Stage1.HIV: Estimation of the community-specific 3yr cumulative HIV incidence
# primary analysis: censors at outmigration/death 
Stage1.HIV <- function(data, subgroup, SL.library=NULL, 
                       outmigrate='Censor', death='Censor'){
  
  # get the W
  W<- get.X(data=data, analysis='HIV')
  
  # Specifying the variables of interest
  E <- c( 'U', 'region_name', 'community_name', 'id', 'pair', 'hiv_prev_0', 
          'mc_cover_0', 'A')
  
  L <- c('outmigrate_3','dead_3', 'oic_hiv_3', 'hiv_3',
         'hivtest_date_0', 'hivtest_date_3')
  
  # Specifying the subgroups if interest
  subgrp <- get.subgroup(data=data, subgroup=subgroup)
 
  OIC <- data$oic_final & subgrp & 
    data$resident_0 & data$stable_0 & data$adult_0 
  
  #-------------------------------------------------------------------
  
  O <- data[OIC, c(E,L)]
  W<- W[OIC, ]
  
  n<- nrow(O)
  
  # time-updated covariates and outcome 
  # assumed to be 0 unless evidence otherwise
  C<-  Delta <- Y <-   rep(0, n)	
  
  # determined to be seroconversion (with full confirmatory testing)
  sc <- O$oic_hiv_3
  
  # non-missing FUY3 status means that completed full algorithm
  Delta[ which( !is.na(sc)) ] <- 1
  
  # primary analysis: all confirmed sercoversions
  Y[ which(sc) ] <- 1
  
  
  # move by FUY3
  move <- which(O$outmigrate_3)
  # dead by FUY3
  dead <- which(O$dead_3)
  
  # if censoring at outmigration (primary)
  if(outmigrate=='Censor'){
    C[move]<- 1
  }
  
  # primary: censoring at death 
  if(death=='Censor'){
    C[dead] <- 1 
  }	else{
    # equivalent to defining the outcome as HIV-free survival
    Delta[dead]<- 1
    Y[dead]<- 1
  }
  
  # Testing time
  time <- as.numeric(O$hivtest_date_3 - O$hivtest_date_0)
  
  
  # individual-level observed data
  O <- data.frame(O[,E], W, C, Delta, Y, time)
  W<- colnames(W)
  
  #----------------------------------------	#----------------------------------------
  # - For each community in turn, we want the counterfactual HIV incidence 
  #		under a hypothetical intervention to prevent right-censoring & 
  #   ensure knowledge of HIV status
  #		E[ Y(c=0, delta=1) ] = P[ Y(c=0, delta=1)= 1]
  # - Primary identifiability:  Y(c, delta) indpt C, Delta 
  # - Secondary identifiability: Y(c, delta) indpt C, Delta | W
  # also return the size of the final OIC
  #----------------------------------------	#----------------------------------------
  
  clusters <- unique(O$id)
  nClust <- length(clusters)
  
  data.clust <- data.frame( matrix(NA, nrow=nClust, ncol=(length(E)+8)  ))
  for(j in 1:nClust){
    
    # now on the OIC
    these.units <- clusters[j]==O$id
    OC <- O[these.units ,]
    
    # the sample size of interest is number in OIC
    oic.size <- nrow(OC)
    
    measured <- OC$C==0 & OC$Delta==1
    sc <- measured &  OC$Y==1
    non.sc <- measured & OC$Y==0
    
    n.outcome <- sum(sc)
    n.measure <- sum(measured)
    
    # sensitivity: person time at risk
    PT <- rep(NA, oic.size)
    PT[sc]<- 0.5*(OC$time[sc])
    PT[non.sc] <- OC$time[non.sc]
    PY.tot<- sum(PT[sc | non.sc])/365
    # incidence eate
    IR<- n.outcome/PY.tot
    
    time.ave <- mean(OC[measured,'time'], na.rm=T)/365
    
    
    # primary outcome: E(Y | C=0, Delta=1) = P(Y =1 | C=0, Delta=1)
    inc.U <- mean(OC[ OC$C==0 & OC$Delta==1, 'Y'])
    # note equivalent to IPTW= mean(  (!OC$C & OC$Delta & OC$Y)/mean(!OC$C & OC$Delta) )
    
    # secondary outcome: E[ E(Y| C=0, Delta=1, W) ]
    Qform<- 'Q.kplus1~ 1'
    OC.sub <- OC[, c(W,'C', 'Delta', 'Y')]
    OC.sub$C <- BinaryToCensoring(is.censored=OC.sub$C)
    inc.A <- ltmle(data=OC.sub, Anodes='Delta', Cnodes='C', 
                     Qform=c(Y=Qform),
                     SL.library=SL.library,
                     Ynodes='Y', abar=1, estimate.time=F, stratify=T)$estimates['tmle']
    
    
    OC$region_name <- as.character(OC$region_name)
    OC$community_name <- as.character(OC$community_name)
    
    # create cluster-level dataframe 
    data.clust[j,] <- c(OC[1,E], n.outcome, n.measure, inc.U, inc.A,
                        oic.size, time.ave, PY.tot, IR)
    print(j)
  }
  colnames(data.clust)<- c(E, 'n.outcome', 'n.measure', 'inc.U', 'inc.A', 
                           'oic.size', 'time.ave', 'PY.tot', 'IR')
  
  list(W=W, data.clust= data.clust)
}



reformat.stage2<- function(data, outcome, indv){
  colnames(data)[grep(outcome, colnames(data) )] <- 'Y'
  colnames(data)[grep(indv, colnames(data))] <- 'nIndv_Y'
  data
}




get.file.name.hiv <- function(subgroup, SL.library, date=NULL){
  
  if(is.null(date)){
    date <- 	format(Sys.time(), "%d%b%Y")
  }
  
  file.name= paste( 'Primary',
                    paste('Subgrp', subgroup, sep=''),
                    paste('v', date,  sep=''), sep="_")
  file.name
}

make.pretty.hiv.min<- function(subgroup='All',
                               date=NULL){
  
  
  file.name= get.file.name.hiv(subgroup, date=date)
  load(paste(file.name, 'Rdata', sep='.') )
  
  # rename subgroups
  if(subgroup=='EU'){
    subgroup <- 'Uganda-East'
  }else if(subgroup=='SWU'){
    subgroup <- 'Uganda-West'
  }else if(subgroup=='Male'){
    subgroup<- 'Men'
  } else if(subgroup=='Female'){
    subgroup <- 'Women'
  } else if(subgroup=='Young'){
    subgroup <- 'Youth (15-24 yr)'
  } else if(subgroup=='Old'){
    subgroup <- 'Aged>24yr'
  } else if(subgroup=='NonMobile'){
    subgroup <- 'Non-mobile'
  } else if(subgroup=='UncircMen'){
    subgroup <- 'Uncircumcised men'
  }
  
  out.csv <- data.frame(matrix(NA, nrow=3, ncol=8))
  colnames(out.csv) <- c('', '', 'pt','CI.lo', 'CI.hi', 'pval', 'N.analysis', 'N.seroc')
  out.csv[,1] <- subgroup
  out.csv[,2] <- c('Intervention', 'Control', 'RR')
  pri <- 'Unadj-Incidence; Adj-Effect'
  out.csv[1, c('pt', 'CI.lo', 'CI.hi')] <- out$est[pri, c('Txt.est', 'Txt.CI.lo', 'Txt.CI.hi')]
  out.csv[2, c('pt', 'CI.lo', 'CI.hi')] <- out$est[pri, c('Con.est', 'Con.CI.lo', 'Con.CI.hi')]
  out.csv[3, c('pt', 'CI.lo', 'CI.hi', 'pval')] <- out$est[pri, c('Effect.est', 'Effect.CI.lo', 'Effect.CI.hi', 'pval')]
  
  txt <- out$data[out$data$A==1,]
  con <- out$data[out$data$A==0,]
  out.csv[1, c('N.analysis', 'N.seroc')] <- c(sum(txt$n.measure), sum(txt$n.outcome))
  out.csv[2, c('N.analysis', 'N.seroc')] <- c(sum(con$n.measure), sum(con$n.outcome))
  out.csv[3, c('N.analysis', 'N.seroc')] <- c(sum(out$data$n.measure), sum(out$data$n.outcome))
  
  out.csv
  
}