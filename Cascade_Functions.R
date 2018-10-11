#==========================================
# R code for stage I estimation in HIV care cascade
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Study Statistician for SEARCH
#==========================================================

Run.Cascade<- function(settings, SL.library='glm', clust.adj){
  
  data.input <- preprocess.cascade()
  
  # Stage 1
  stage1  <- Stage1.Cascade.Xsect(data.input, settings=settings,
                                  SL.library=SL.library)
  covariates <- stage1$W
  data <- stage1$data.clust
  raw2 <- stage1$RAW2
  primary <- stage1$primary
  secondary <- stage1$secondary
  
  est <- arm.pri1 <- arm.pri0 <- raw2.sec1 <- raw2.sec0 <- NULL
  
  
  # primary analysis weights wrt estimated number of HIV+
  alpha.pri <- get.weights(primary, weighting='indv')$alpha
  alpha.sec <- get.weights(secondary, weighting='indv')$alpha
  # weighting prevalence wrt population size
  temp <- data[,c('id', 'N.pop')]
  colnames(temp) <- c('id', 'nIndv')
  alpha.prev <- get.weights(temp, weighting='indv')$alpha
  # calculate arm-specific cascade
  arm.pri <- get.subgroup.est(data=primary, alpha= alpha.pri, alpha.prev=alpha.prev)
  arm.sec <- get.subgroup.est(data=secondary, alpha=alpha.sec, alpha.prev=alpha.prev)
  
  if( var(data$A)!=0 ){    # then varability in the exposure
    
    # STAGE 2 ANALYSIS:
    # First with the primary Yc outcome 
    data.temp <- reformat.stage2(primary)
    Adj.supp_Adj.effect <- Stage2(weighting='indv', data.input=data.temp, 
                                  outcome='Y', clust.adj=clust.adj,
                                  do.data.adapt=T)
    
    # secondary yc
    data.temp <- reformat.stage2(secondary)
    Unadj.supp_Unadj.effect <- Stage2(weighting='indv', data.input=data.temp, 
                                      outcome='Y',  do.data.adapt=F)
    est <- data.frame(rbind(Adj.supp_Adj.effect=Adj.supp_Adj.effect, 
                            Unadj.supp_Unadj.effect=Unadj.supp_Unadj.effect))
    txt<- primary$A==1
    con<- primary$A==0
    arm.pri1 <- get.subgroup.est(data=primary[txt,], alpha= alpha.pri[txt],  
                                 alpha.prev=alpha.prev[txt])
    arm.pri0 <- get.subgroup.est(data=primary[con,], alpha= alpha.pri[con], 
                                 alpha.prev=alpha.prev[con])
        # unadjusted output
    raw2.sec1 <- get.supp.table(raw2=raw2[raw2$A==1,], pri=arm.pri1)
    raw2.sec0 <- get.supp.table(raw2=raw2[raw2$A==0,], pri=arm.pri0)
  
  }
    

  rownames(data) <- rownames(primary) <- rownames(secondary) <- data$community_name
  
  data<- subset(data, select=-c(U,community_name, id))
  data[data$region_name=='Western Uganda','region_name'] <- 3
  data[data$region_name=='Eastern Uganda', 'region_name'] <- 2
  data[data$region_name=='Kenya', 'region_name'] <- 1
  colnames(data)[1]<- 'region'
  
  data$region <- as.numeric(data$region)
  
  primary<- subset(primary, select=-c(U,region_name, community_name, id))
  secondary<- subset(secondary, select=-c(U,region_name, community_name, id))
  
  list(covariates=covariates, data=data, primary=primary, secondary=secondary, est=est,
       arm.pri = arm.pri, arm.sec = arm.sec,
       arm.pri1 = arm.pri1, raw2.sec1=raw2.sec1,
       arm.pri0 = arm.pri0, raw2.sec0=raw2.sec0)
  
}



preprocess.cascade <- function(){
	
	load("outputs-withIntOnly.RData")
	print( dim(outputs) )
	data.input <- outputs

	# Exclude anyone that was flagged as an SEARCH-id related error
	data.input <- subset(data.input, !(data_flag | dead_0 | move_0) )
	
	# outputs-withIntOnly.Rdata does NOT have community_number
	#		add now
	id<-  rep(NA, nrow(data.input ))
	
	comm <- unique(data.input$community_name)
	for(j in 1:32){
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


# Stage1.Cascade.Xsect: get the community-specific estimates of VL suppression
#	in open cohort of HIV+ at t
Stage1.Cascade.Xsect <- function(data.input, settings,
                                 SL.library){
  
  # cluster-level covariates
  E <- c( 'U', 'region_name', 'community_name', 'id', 'pair', 'A') 
  
  # define target population for the analysis
  restrict<- get.pop.cascade(data=data.input, settings=settings)
  
  data.input <- subset(data.input, restrict)
  print(dim(data.input))
  
  clusters <- unique(data.input$id)
  nClust <- length(clusters)
  
  data.clust <- data.frame(matrix(NA, nrow=nClust, ncol=(length(E) + 12  )))
  colnames(data.clust) <- c(E, "N.pop", "N.chc", "N.trk", "N.know.HIV",
                            "N.Delta", "N.know.pos", "N.est.pos" ,"N.know.pDx",
                            "N.know.eART", "N.know.supp", "N.TstVL", "N.est.supp")
  
  primary <- data.frame(matrix( NA, nrow=nClust, ncol=(length(E) + 19)))
  colnames(primary)<- c(E, 'pSupp_0', 'pAdol_0', 
                        "prev_pt","prev_CI.lo", "prev_CI.hi",
                        "pDx.pos_pt", "pDx.pos_CI.lo", "pDx.pos_CI.hi",
                        "eART.pDx_pt", "eART.pDx_CI.lo", "eART.pDx_CI.hi",
                        "supp.eART_pt", "supp.eART_CI.lo", "supp.eART_CI.hi",
                        "supp.pos_pt", "supp.pos_CI.lo","supp.pos_CI.hi",
                        "Product", "nIndv")
  secondary <- primary 
  
  # unadjusted numbers for supplementary table
  RAW2 <- data.frame(matrix(NA, nrow=nClust, ncol=length(E)+10))
  
  for(j in 1:nClust){
    
    these.units <- data.input$id==clusters[j]
    OC <- data.input[these.units ,]
    
    baseline.pred <-  get.X(data=OC, analysis='Cascade', 
                            adj.full=settings$adj.full) 
    # print(names(baseline.pred))
    out<- do.serial.analysis(OC= OC, baseline.pred=baseline.pred,
                             settings=settings,
                             SL.library= SL.library)
    
    # aggregate relevant covariates
    OC$region_name <- as.character(OC$region_name)
    OC$community_name <- as.character(OC$community_name)
    data.clust[j,] <- c(OC[1,E], out$RAW) 
    RAW2[j,] <- c(OC[1,E], out$RAW2)
    
    # create cluster-level dataframe 
    pSupp_0 <- data.frame(pSupp_0=
                            mean(OC[OC$hiv_0==1, 'supp_0'], na.rm=T))
    pAdol_0 <- data.frame(pAdol_0=
                            mean(OC[OC$hiv_0==1, 'age_0']<=24, na.rm=T))
    
    # create cluster-level dataframe 
    
    primary[j,] <- c(OC[1,E],  pSupp_0, pAdol_0, out$Cascade, 
                     nIndv=out$RAW$N.est.pos)
    secondary[j,] <- c(OC[1,E],pSupp_0, pAdol_0, out$Cascade.sec,
                       # number estimated to be HIVpositve
                       nIndv=round(out$RAW$N.pop*out$Cascade.sec$prev_pt))
    print(j)                
  }
  colnames(RAW2) <- c(E, colnames(out$RAW2))
  
  
  W <- names(baseline.pred)
  list(W=W, data.clust=data.clust, RAW2=RAW2,
       primary=primary,
       secondary=secondary )
}



get.pop.cascade <- function(data, settings){
	
	n<- nrow(data)
	
	region <- get.subgroup(data=data, subgroup=settings$region)

	# Specifying the subgroups if interest
	this.subgroup <- get.subgroup(data=data, subgroup=settings$subgroup, 
	                              time=settings$time)
	
	# have to be 15+, cannot have died or outmigrated at t 
	time <- settings$time
	if(time==0){
		adult <- data$age_0>14
		alive <-  rep(T, n)
		move <- rep(F,n)
	} else if(time==1){
		adult <- data$age_0>13
		alive <- !data$dead_1
		move <- data$move_1
	} else if (time==2){
		adult <- data$age_0>12
		alive <- !data$dead_2
		move <- data$move_2
	} else if (time==3){
		adult <-  rep(T, n) 
		dead <- move <- rep(F, n)
		dead[which(data$dead_3)] <- T
		alive <- !dead
		move[which(data$outmigrate_3)] <-T
	}
	
	if(time==1 | time==2){
		arm <- data$A==1
	} else{
		arm <- rep(T, n)
	}
	
	
	if(time==0){
	  resident <- data$resident_0
	} else if (time==1){ # Y1
	  resident <- data$resident_recensus_1   
	} else if(time==2){
	  resident <- data$resident_recensus_2
	} else if(time==3){ 
	  resident <- data$resident_3
	}
	
	restrict<- region &  this.subgroup & adult & alive & !move & arm & resident 
	restrict
}


#*-------


# do.serial.analysis - code to run primary and secondary analyses for 
#	serial cross-sectional analysis of prevalent HIV+ 
do.serial.analysis <- function(OC, baseline.pred, settings, 
                               SL.library = NULL, verbose = F) {


  # used in deterministic Q adjustment sets
	data_0 <- preprocess.serial(data = OC, time=0, settings=settings)
	colnames(data_0) <- paste(colnames(data_0), "0", sep = "_")

  if(settings$time==2){
  	data_1 <- preprocess.serial(data = OC, time=1, settings=settings)
	  colnames(data_1) <- paste(colnames(data_1), "1", sep = "_")
  } else{
    data_1 <-  NULL
  }  
	
	time <- settings$time
	data <- preprocess.serial(data = OC, time=time, settings=settings)


	#**************** Prevalence of HIV at time t: P(Y*=1)
	Prob.HIVpos <- get.prevalence(baseline.pred=baseline.pred, 
	                              data_0=data_0, data_1=data_1, 
	                              data=data, time=time, 
	                              SL.library=SL.library, verbose=verbose)
	
	
	#****************   Previous diagnosis at time t: P(pDx=1, Y*=1) 
	Prob.pDx <- get.pDx(baseline.pred=baseline.pred, 
	                    data_0=data_0, data_1=data_1,
	                    data=data, time=time, 
	                    SL.library=SL.library, verbose=verbose)
  
  

	
	#**************************** Ever ART use at time t: P(eART=1, pDx=1, Y*=1) 
	Prob.eART <- get.eART( baseline.pred=baseline.pred, 
	                      data_0=data_0, data_1=data_1,
	                      data=data, time=time,
	                      SL.library=SL.library, verbose=verbose)
	
	#****************Suppression at t:  P(Supp*=1, eART=1, pDx=1, Y*=1) 
	Prob.joint <- get.joint(baseline.pred=baseline.pred,
	                        data_0=data_0, data_1=data_1,
	                        data=data, time=time,
	                        SL.library=SL.library, verbose=verbose) 
	

	
	#************************** Compiling the results
	# Number in the population of interest
	N.pop <- nrow(data)
	RAW <- data.frame(
	  # Number in the population of interest
	  N.pop=N.pop,
	  # Number seen at CHC
	  N.chc = sum(data$chc),
	  N.trk = sum(data$tr),
	  N.know.HIV = sum(data$TstHIV),
	  # number seen with known status at CHC/track
	  N.Delta = sum(Prob.HIVpos$sum.A),
	  # Number known to be HIV+
	  N.know.pos = sum(data$HIVpos),
	 	# Number estimated to be HIV+ = (Estimated HIV prevalence) x (Population size) 
	  N.est.pos = round(Prob.HIVpos$pri$e$pt*N.pop, 0),
	  # Number with prior diagnosis in the population
	  N.know.pDx = round(Prob.pDx$pri$e$pt*N.pop, 0),
	  # Number with eART in the population
	  N.know.eART = round(Prob.eART$pri$e$pt*N.pop, 0),
	  N.know.supp = sum(data$Supp),
	  # Number with relevant VL measure
	  N.TstVL=sum(Prob.joint$sum.A),
	  # Number estimated to be suppressed = (Estimated Suppression ) x (Population size)
	  N.est.supp =round(Prob.joint$pri$e$pt * N.pop, 0)
	)


	Cascade <- do.ratios(Prob.HIVpos = Prob.HIVpos$pri, 
	                     Prob.pDx = Prob.pDx$pri, 
	                     Prob.eART = Prob.eART$pri, 
	                     Prob.joint = Prob.joint$pri, 
	                     primary = T)

	Cascade.sec <- do.ratios(Prob.HIVpos = Prob.HIVpos$sec, 
	                         Prob.pDx = Prob.pDx$sec, 
	                         Prob.eART = Prob.eART$sec, 
	                         primary = F, data = data, verbose = verbose)

	
	# should correspond to unadjusted estimates 
	RAW2 <- data.frame(
	  # all known HIV+ regardless if seen at FUY3
	  allpos =  sum(data$HIVpos),
	  # all known HIV+ with missed VL - regardless if seen
	  allpos.noTstVL = sum(data$HIVpos*!data$TstVL),
	  #
	  # all folling conditions on being seen at CHC/tracking (Delta=1)
	  # prevalence = A:delta, Y:HIVpos 
	  HIVpos  = sum(data$HIVpos*data$Delta),
  	# pDx: A=:Delta, Y=pDx
    pDx =sum(data$pDx*data$Delta),
	  # eART: A=Delta, Y=eART
	  eART = sum(data$eART*data$Delta),
  	# with VL measured: A=Delta*TstVL Y = eART
	  eART.TstVL = sum(data$eART*data$Delta*data$TstVL),
	  # supp with ART: A=Delta*TstVL Y = supp
  	supp.eART.TstVL = sum(data$Supp*data$eART*data$Delta*data$TstVL),
	  # supp HIV+: A=Delta*TstVL Y = supp
	  supp.TstVL = sum(data$Supp*data$Delta*data$TstVL),
    # HIV+ with measured VL
    HIVpos.TstVL =sum(data$HIVpos*data$Delta*data$TstVL),
    # HIV+ with miss VL
    HIVpos.NoTstVL = sum(data$HIVpos*data$Delta*!data$TstVL) )

	

	
	list(RAW = RAW, RAW2=RAW2, Cascade = Cascade, 
		Cascade.sec = Cascade.sec)

}



preprocess.serial <- function(data, time, settings){
  
  if(time==0 | time==3){
    # if using parallel data in intervention and control
    HIV.variable <- 'hiv'
    pDx.variable <- 'pdx_vl'
    eART.variable <- 'eart_vl'
  } else{
    # if, instead, using interim data in the intervention arm
    HIV.variable	<- 'all_hiv'
    pDx.variable <- 'hiv_preCHC'
    eART.variable <- 'i_eart_vl'  
  }
  
  
  n <- nrow(data)
  
  pDx <- eART <- Delta  <- HIVpos <- TstVL <- Supp <- rep(0, n)
  
  # no evidence of prior Dx or ART == fail
  pDx[ which(data[, paste(pDx.variable, time, sep='_')] ==1) ] <-1 
  eART[ which(data[, paste(eART.variable, time, sep='_')] ==1)] <-1 
  # if on ART then previously dx-ed
  pDx[ eART==1] <- 1
  
  chc <- as.numeric(data[, paste('chc', time, sep='_')]	)
  tr <- as.numeric(data[, paste('tr', time, sep='_')])

  # did we actually know your HIV status?
  hiv.temp <-  data[, paste(HIV.variable, time, sep='_')]	
  TstHIV <- as.numeric( !is.na(hiv.temp) )
  # Delta - require that we saw you at t & had a known status
  Delta[ which( (chc==1 | tr==1)  & TstHIV==1 )] <-1
  
  #  HIV status as =1 if HIV+ and 0 otherwise 
  HIVpos[ which(hiv.temp) ] <- 1
  
  # Suppression - VL.variable='supp' 
  supp.temp <- data[, paste('supp', time, sep='_')]
  # SET SUPP=NA if not HIVpos
  supp.temp[ which( !is.na(supp.temp) & HIVpos==0) ] <- NA
  TstVL[ which(!is.na(supp.temp) ) ] <- 1
  
  #  suppression as =1 if suppressed and 0 otherwise
  Supp[ which(supp.temp) ] <- 1
  
  data.frame(cbind(pDx, eART, chc, tr, TstHIV, Delta, HIVpos, TstVL, Supp))
}


#**************** Prevalence of HIV at time t: P(Y*=1)
get.prevalence <- function(baseline.pred, data_0, data_1,  data,
                           time, SL.library, verbose=F){
  
  # outcome as observed HIV status 
  Y <- data$HIVpos 
  # intervention variable: seen at CHC/tracking with HIV status known
  A <- data$Delta 
  
  # specify adjustment set other than baseline covariates
  adj <- get.adjustment(data = data,
                        baseline.pred = baseline.pred,
                        data_0 = data_0, data_1=data_1,
                        time = time)
  
  # using deterministic knowledge about known HIV status at prior time point 
  Prob.HIVpos <- call.ltmle(pred.A = adj$adj, 
                            A = A, Y = Y, SL.library = SL.library, 
                            deterministicQ = adj$detQ)
  if(verbose)print(Prob.HIVpos$est)
  
  #* Secondary analysis 
  Prob.HIVpos.sec <- do.secondary(A = A, Y = Y, verbose = verbose)
  
  list(sum.A=sum(A), pri=Prob.HIVpos, sec=Prob.HIVpos.sec)
  
}


#****************   Previous diagnosis at time t: P(pDx=1, Y*=1) 
get.pDx <- function(baseline.pred, data_0, data_1, data,
                    time, SL.library, verbose=F){
  
  Y <- data$pDx

  # assumes complete capture
  Prob.pDx <- call.ltmle(pred.A = NULL, 
                         A = rep(1, nrow(data)), Y = Y, 
                         SL.library = NULL)
  # c(Prob.pDx$est$pt, sum(Y)/length(Y))
  if(verbose)print(Prob.pDx$est)
  
  #* Secondary analysis 
  Prob.pDx.sec <- do.secondary(A = data$Delta,
                               Y = Y, verbose = verbose)
  list(pri=Prob.pDx, sec=Prob.pDx.sec)
}




#**************************** Ever ART use at time t: P(eART=1, pDx=1, Y*=1) 

get.eART <-  function(baseline.pred, data_0, data_1,
                      data, time, SL.library, verbose=F){
  Y <- data$eART
 
  Prob.eART <- call.ltmle(pred.A = NULL,
                            A = rep(1, nrow(data)), Y = Y, 
                            SL.library = NULL)
  # c(Prob.eART$est$pt, sum(Y)/length(Y))
  if(verbose)print(Prob.eART$est)
  
  #* Secondary analysis  
  Prob.eART.sec <- do.secondary(A = data$Delta, 
                                Y = Y, verbose = verbose)

  list(pri=Prob.eART, sec=Prob.eART.sec)
}


#****************Suppression at t:  P(Supp*=1, eART=1, pDx=1, Y*=1) 
get.joint <- function(baseline.pred, data_0, data_1,
                      data, time, SL.library, verbose=F){  
  
  # last QC that if HIV+&supp, must have started ART
  if((sum(data$eART==0 & data$Supp==1)>0)){
    Y <- data$eART*data$Supp
  }else{
    Y <- data$Supp
  }
  A <- data$TstVL 
  
  # get adjustment set using deterministic knowledge (primary)
  pred.A <- get.adjustment.joint(data = data,
                                 baseline.pred = baseline.pred,
                                 data_0 = data_0, data_1=data_1,
                                 time = time)
  
  Prob.joint <- call.ltmle(pred.A = pred.A, 
                             A = A, Y = Y, SL.library = SL.library, 
                             deterministicQ = deterministicQ_NO)
  if(verbose)print(Prob.joint$est)
  
  list(sum.A=sum(A), pri=Prob.joint)
}  



# get.adjustment: function to get the adjustment set for prevalence
get.adjustment<- function(data, baseline.pred, 
                          data_0, data_1,
                          time){
  
  # always adjusting for baseline testing
  TstHIV_0 <- data_0$TstHIV_0
  
  if(time==0){
    adj <- baseline.pred
    detQ<- NULL
  } else { 
    if(time==1 | time==3){
      # at Y1 or Y3 (not using interim data to be parallel across arms)
      adj <- data.frame(baseline.pred, TstHIV_0, 
                        detQ.variable=data_0[, paste('HIVpos', 0, sep='_')])
    }else{
      # if time=2 adjust for baseline and time-1
      adj <- data.frame(baseline.pred, TstHIV_0, 
                          detQ.variable=data_1[, paste('HIVpos', 1, sep='_')])
    }
    detQ <- deterministicQ_YES
  }
  list(adj=adj, detQ=detQ)
}



get.adjustment.joint <- function(data, baseline.pred, data_0, 
                                 data_1, time){
  
  TstHIV_0 <- data_0$TstHIV_0
  
  if(time==0){
    adj <- baseline.pred
  } else if( time==1 | time==3 ){
    adj <- data.frame(baseline.pred, TstHIV_0,
                      prior=data_0[, paste('Supp', 0, sep='_')])
  }else{
    # if time=2 adjust for baseline and time=1
    adj <- data.frame(baseline.pred, TstHIV_0,
                      prior=data_1[, paste('Supp', 1, sep='_')])
    
  } 
  
  adj <- data.frame(adj, detQ.variable=data$eART)
  
  
  adj
  
}



# call.ltmle: function to call the ltmle 
call.ltmle<- function(pred.A=NULL, A,  Y, 
                      SL.library=NULL,
                      deterministicQ=NULL, 
                      observation.weights=NULL, 
                      id=NULL, 
                      verbose=F){
  
  
  # create temporary data frame
  if( is.null(pred.A) ){
    data.temp<- data.frame(A, Y)
  } else{
    data.temp<- data.frame(pred.A, A, Y)
  }

  est.temp<- ltmle(data=data.temp, Anodes='A', Lnodes=NULL, Ynodes='Y',
                   abar=1,
                   stratify = T, SL.library=SL.library,
                   estimate.time=F,
                   variance.method='ic', 
                   deterministic.Q.function= deterministicQ,
                   observation.weights=observation.weights, id=id)
  
  
  if(verbose){
    print(est.temp$fit$g)
    print(est.temp$fit$Q)
  }				
  IC<- est.temp$IC$tmle
  
  est<- data.frame(pt=est.temp$estimate["tmle"], 
                   CI.lo=summary(est.temp)$treatment$CI[1], 
                   CI.hi=summary(est.temp)$treatment$CI[2] 	)
  
  list(est=est,IC=IC)
}

# simple function do the secondary analysis
do.secondary <- function(A, Y, verbose=F){
  N.measured <- sum(A)
  N.obs.outcome <- sum(A & Y)
  
  out<-  call.ltmle(pred.A= NULL, A=A, Y=Y, SL.library=NULL)
  if(verbose){
    print(c(out$est$pt,N.obs.outcome/ N.measured))
  }
  out
}

#*====
# USE THE ABOVE FOUR PROBABILITIES TO CALCULATE CASCADE COVERAGE 
# 	USE THE DELTA.METHOD FOR INFERENCE IN PRIMARY ANALYSIS
#   See get.var.bayes.bayes function
do.ratios <- function(Prob.HIVpos, Prob.pDx, Prob.eART, Prob.joint,
                      primary=T, data, verbose=F){
  
  prev<- Prob.HIVpos$est
  colnames(prev) <- paste('prev', colnames(prev), sep='_')
  
  # calculate
  pDx.pos <- get.var.bayes(mu1=Prob.pDx$est$pt, 
                           IC1=Prob.pDx$IC, 
                           mu0=Prob.HIVpos$est$pt,
                           IC0=Prob.HIVpos$IC)$est
  colnames(pDx.pos) <- paste('pDx.pos', colnames(pDx.pos), sep='_')
  
  eART.pDx <- 	get.var.bayes(mu1=Prob.eART$est$pt,
                             IC1=Prob.eART$IC, 
                             mu0=Prob.pDx$est$pt,
                             IC0=Prob.pDx$IC)$est
  colnames(eART.pDx) <- paste('eART.pDx', colnames(eART.pDx), sep='_')
  
  if(primary){
    supp.eART <- get.var.bayes(mu1=Prob.joint$est$pt,
                               IC1=Prob.joint$IC,
                               mu0=Prob.eART$est$pt, 
                               IC0=Prob.eART$IC)$est
    
    supp.pos <- get.var.bayes(mu1=Prob.joint$est$pt,
                              IC1=Prob.joint$IC, 
                              mu0=Prob.HIVpos$est$pt,
                              IC0=Prob.HIVpos$IC)$est 
  }else{
    # identifiability requires only using VL at CHC/tracking
    data.sec <- data
    
    data.sec$TstVL <- data.sec$TstVL*data.sec$Delta   
    data.sec$Supp <- data.sec$Supp*data.sec$TstVL
    
    on.ART<- data.sec$eART==1
    supp.eART <- do.secondary(A=data.sec$TstVL[on.ART],
                              Y=data.sec$Supp[on.ART],
                              verbose=verbose)$est
    
    supp.pos <- do.secondary(A=data.sec$TstVL, 
                             Y=data.sec$Supp, 
                             verbose=verbose)$est
  }
  colnames(supp.eART) <- paste('supp.eART', colnames(supp.eART), sep='_')
  colnames(supp.pos) <- paste('supp.pos', colnames(supp.pos), sep='_')
  
  Product <- as.numeric(pDx.pos[1])*
    as.numeric(eART.pDx[1])*
    as.numeric(supp.eART[1])
  data.frame(prev, pDx.pos, eART.pDx, supp.eART, supp.pos, Product)
}


#*====

#* STAGE 2

reformat.stage2<- function(data){
  colnames(data)[grep('supp.pos_pt', colnames(data) )] <- 'Y'
  colnames(data)[grep('nIndv', colnames(data))] <- 'nIndv_Y'
  data
}

get.subgroup.est<- function(data, alpha, alpha.prev){	
  
  A <- rep(1, nrow(data))
  # weight wrt number of indv in population
  prev <- call.ltmle(A=A, Y=data$prev_pt, observation.weights=alpha.prev)$est
  # weight following wrt number of estimated HIV+
  pDx.pos <- call.ltmle(A=A, Y=data$pDx.pos_pt, observation.weights=alpha)$est
  eART.pDx <- call.ltmle(A=A, Y=data$eART.pDx_pt,  observation.weights=alpha)$est
  supp.eART <- call.ltmle(A=A, Y=data$supp.eART_pt,  observation.weights=alpha)$est
  supp.pos<- call.ltmle(A=A, Y=data$supp.pos_pt, observation.weights=alpha)$est
  Product<- call.ltmle(A=A, Y=data$Product, observation.weights=alpha)$est
  
  data.frame(rbind(prev=prev, pDx.pos=pDx.pos, eART.pDx= eART.pDx, 
                   supp.eART= supp.eART,  supp.pos= supp.pos, Product=Product))
  
}


get.file.name.cascade <- function(region='All', subgroup='All', 
                                   time=3, 
                                  date=NULL){
	
	if(is.null(date)){
		date <- 	format(Sys.time(), "%d%b%Y")
	}
  adj.full <- subgroup!='Young'
	file.name= paste( 'Cascade',
	                  region,
		subgroup,
		paste('Year', time, sep=''),
		paste('FullAdjust', adj.full, sep=''),
		paste('v', date,  sep=''), sep="_")
	file.name
}

#*====



# # alternative calculation for the unadjusted/secondary estimates
supp.table.helper <- function(num, den){
  num <- sum(num);
  den <- sum(den)
  paste( round(num/den*100), '% (', paste(num, den, sep='/'), ')', sep='')
}
supp.table.helper2 <- function(pri, this){
  paste( round(pri[this, 'pt']*100), '% (',
         round(pri[this, 'CI.lo']*100), ',',
         round(pri[this, 'CI.hi']*100), ')', 
         sep='')
         
}
get.supp.table <- function(raw2, pri){
  out<- data.frame(rbind(
    pDx.pos= c( supp.table.helper(num=raw2$pDx, den=raw2$HIVpos),
                supp.table.helper2(pri, this='pDx.pos') ),
    eART.pDx= c(supp.table.helper(num=raw2$eART,  den=raw2$pDx),
                supp.table.helper2(pri, this='eART.pDx') ),
    supp.eART= c(supp.table.helper(num=raw2$supp.eART.TstVL,
                                   den=raw2$eART.TstVL),
                 supp.table.helper2(pri, this='supp.eART') ),
    supp.pos = c(supp.table.helper(num=raw2$supp.TstVL,  
                                   den=raw2$HIVpos.TstVL),
                 supp.table.helper2(pri, this='supp.pos') ),
    pos.noTstVL = c( supp.table.helper(num=raw2$HIVpos.NoTstVL,  
                                       den=raw2$HIVpos),
                     NA),
    allpos.noTstVL = c( supp.table.helper(num=raw2$allpos.noTstVL,  
                                       den=raw2$allpos),
                     NA)
    ))
   
  colnames(out) <- c('U: % (num/den)', 'A: % (CI in %)' )
  out
}


make.pretty.min <- function(out, region='All', subgroup='All', 
                            time=3, date='final'){
  file.name<- get.file.name.cascade(region=region,
                             subgroup=subgroup, 
                             time=time, 
                             date=date)
  load(paste(file.name, 'Rdata', sep='.'))
  if(subgroup=='All'){
    subgroup<- region
  } else{
    subgroup <- subgroup
  }
  effect <- out$est['Adj.supp_Adj.effect',]
  pool <- out$arm.pri
  arm1 <- out$arm.pri1
  arm0 <- out$arm.pri0
  
  out.csv <- data.frame(matrix(NA, nrow=3, ncol=9))
  colnames(out.csv) <- c('', '', '','pt','CI.lo', 'CI.hi', 'pval', 'N.analysis', 'N.est.HIV+')
  out.csv[,1] <-  subgroup
  out.csv[,3] <- c('Intervention', 'Control', NA)
  
  basic<- c('pt', 'CI.lo', 'CI.hi')
  
  # GET N   
  txt <- out$data[out$data$A==1,]
  con <- out$data[out$data$A==0,]
  out.csv[1, c('N.analysis', 'N.est.HIV+')] <- c(sum(txt$N.pop), 
                                                 round( sum(txt$N.pop)*arm1['prev','pt']))
  out.csv[2, c('N.analysis', 'N.est.HIV+')] <- c(sum(con$N.pop), 
                                                 round( sum(con$N.pop)*arm0['prev','pt']))
  out.csv[3,c('N.analysis', 'N.est.HIV+')]  <- c(sum(out$data$N.pop), NA)
  
  if(settings$time==0){ 
    #IGNORE EFFECT
    out.csv[,2] <- 'Baseline';   out.csv[3,3] <- 'Overall'
    out.csv[1, basic] <- arm1['supp.pos', c('pt', 'CI.lo', 'CI.hi')]
    out.csv[2, basic] <- arm0['supp.pos', c('pt', 'CI.lo', 'CI.hi')]
    out.csv[3, basic] <- pool['supp.pos', c('pt', 'CI.lo', 'CI.hi')]
    
  }else if (settings$time==3){
    out.csv[,2] <- 'Year 3'
    out.csv[1, basic] <- effect[, c('Txt.est', 'Txt.CI.lo', 'Txt.CI.hi')]
    out.csv[2, basic] <- effect[, c('Con.est', 'Con.CI.lo', 'Con.CI.hi')]
    
    out.csv[3, c(basic, 'pval')] <- effect[, c('Effect.est', 'Effect.CI.lo', 'Effect.CI.hi', 
                                               'pval')]
    out.csv[3,3] <- 'RR'
    
    
  } else{
    if(settings$time==1){
      out.csv[,2] <- 'Year 1'
    } else if(settings$time==2){
      out.csv[,2] <- 'Year 2'
    } 
    out.csv[, basic] <- pool['supp.pos', basic]
    out.csv[, c('N.analysis', 'N.est.HIV+')] <- c(sum(out$data$N.pop), 
                                                  round(sum(out$data$N.pop)*pool['prev','pt']))
    out.csv<- out.csv[1,]
  }
  
  
  out.csv
  
}



