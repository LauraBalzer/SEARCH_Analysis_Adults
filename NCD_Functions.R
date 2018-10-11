#==========================================
# R code for stage I estimation in NCDs
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Study Statistician for SEARCH
#==========================================================


Run.NCD <- function(settings, SL.library='glm'){
  
  data.input<- get.full.data.ncd(settings=settings)

  # outcome A: All* with prevalent HT/NCD at FUY3 
  outA  <- Stage2.NCD(data.input=data.input, settings=settings, outcome='A',
                        SL.library=SL.library)
  prev.control <- get.CI.prev.control(ncd= settings$ncd,
                                      data.clust=outA$data.clust, 
                                      weighting='indv') 

  if(settings$time==0){
    # then only getting baseline prevalence estimates
    RETURN <- prev.control
  }else{
    # outcome B: HIV+* with prevalent HT/NCD at FUY3
    outB  <- Stage2.NCD(data.input=data.input, settings=settings, outcome='B',
                        SL.library=SL.library)
    
    prev.control.HIVpos <- get.CI.prev.control(ncd= settings$ncd,
                                               data.clust=outB$data.clust, 
                                               weighting='indv') 
    
    # outcome C: HIV+:Dual control among prevalent HT/NCD at FUY3**  
    outC  <- Stage2.NCD(data.input=data.input, settings=settings, outcome='C',
                        SL.library='glm')
    
    # outcome D: All with prevalent HT/NCD at BL  
    outD  <- Stage2.NCD(data.input=data.input, settings=settings, outcome='D',
                           SL.library=SL.library)
    
    # outcome E: HIV+ and prevalent HT/NCD at BL  
    outE  <- Stage2.NCD(data.input=data.input, settings=settings, outcome='E',
                          SL.library='glm')
    
    # outcome F: HIV+:Dual control among prevalent HT/NCD at BL**   
    outF  <- Stage2.NCD(data.input=data.input, settings=settings, outcome='F',
                          SL.library='glm')

    EST <- rbind(outA$est,
                 outB$est,
                 outC$est,
                 outD$est,
                 outE$est,
                 outF$est
    )
    rownames(EST) <- c(
      'A: All* with prevalent ncd at Y3; Adj.',
      'A: All* with prevalent ncd at Y3; Unadj',
      'B: HIV+* with prevalent ncd at Y3; Adj.',
      'B: HIV+* with prevalent ncd at Y3; Unadj.',
      'C: HIV+:Dual control among prevalent ncd at Y3**; Adj.',
      'C: HIV+:Dual control among prevalent ncd at Y3**; Unadj.',
      'D: All with prevalent ncd at BL; Adj.',
      'D: All with prevalent ncd at BL; Unadj.',
      'E: HIV+ and prevalent ncd at BL; Adj.',
      'E: HIV+ and prevalent ncd at BL; Unadj.',
      'F: HIV+:Dual control among prevalent ncd at BL**; Adj.',
      'F: HIV+:Dual control among prevalent ncd at BL**; Unadj'
    )
    
    
    
    RETURN <- list(EST=EST,  
                   outA=outA$data.clust,
                   outB=outB$data.clust,
                   outC=outC$data.clust,
                   outD=outD$data.clust, 
                   outE=outE$data.clust,
                   outF=outF$data.clust, 
                   W=outA$W, W.red=outC$W, 
                   prev.control=prev.control,
                   prev.control.HIVpos=prev.control.HIVpos)
  }

  RETURN
  

}




#PREPROCESS
get.full.data.ncd <- function(settings){
  
  # load complete dataset
  load("outputs-withIntOnly.RData")
  print( dim(outputs) )
  data.input <- outputs
  
  # Exclude anyone that was flagged as an SEARCH-id related error
  data.input <- subset(data.input, !(data_flag | dead_0 | move_0) )
  
  restrict <- get.restrict.ncd(data.input=data.input, 
                               settings=settings)
  data.input <- data.input[restrict,]

  # outputs-withIntOnly.Rdata does NOT have community_number
  # get the cluster-level adjustment variables
  id <- ncd_control_0 <- dual_control_0 <- chc_cover_0 <- 
    alcohol_prev_0 <- overwt_prev_0 <-  rep(-99, nrow(data.input ))
  
  comm <- unique(data.input$community_name)
  
  for(j in 1:length(comm)){
    these.units <- data.input$community_name==comm[j]
    id[these.units] <- j
    
    OC<- data.input[these.units, ]

    #  BL control among BL prev (unadjusted)
    temp <-  get.pop(data=OC, ncd=settings$ncd, outcome='D')
    BL <-  preprocess.NCD(data=temp, time=0, ncd=settings$ncd, outcome='D')
    ncd_control_0[these.units]  <- sum(BL$delta & BL$Y)/sum(BL$delta) 		
    
    # BL dual-control among BL NCD-HIV+ (unadjusted)
    temp <-  get.pop(data=OC, ncd=settings$ncd, outcome='F')
    BL <- preprocess.NCD(data=temp, time=0, ncd=settings$ncd, outcome='F')
    dual_control_0[these.units]  <- sum(BL$delta & BL$Y)/sum(BL$delta) 
       
    chc_cover_0[these.units] <- mean(OC$chc_0, na.rm=T)
    alcohol_prev_0[these.units] <- mean(OC$alcohol_0, na.rm=T)
    overwt_prev_0[these.units] <- mean( OC$bmi_0>=25, na.rm=T)
  }	
  
  data.input$region_name <- as.character(data.input$region_name)
  data.input$community_name <- as.character(data.input$community_name)
  
  
  data.input <- cbind(data.input, 
                      A= as.numeric(as.logical(data.input$intervention)),
                      id=id,
                      ncd_control_0=ncd_control_0,
                      dual_control_0=dual_control_0,
                      chc_cover_0=chc_cover_0,
                      alcohol_prev_0=alcohol_prev_0,
                      overwt_prev_0=overwt_prev_0
                      )
  
  # transform pairs to be numeric
  data.input$pair <- as.numeric(as.character(data.input$pair))
  
  # if haven't added a dummy variable for unadjusted
  if( sum(grep('U', colnames(data.input)))==0){
    data.input <- cbind(U=1, data.input) 
  }
  
  print( dim(data.input) )
  
  print('***preprocessing done***')
  
  data.input
}


get.restrict.ncd<- function(data.input, settings){
   
  
  # age restriction
  if(settings$time==0){
    age <- data.input$age_0 >29
  } else if(settings$time==3){
    age <- data.input$age_3 >29
  }
  
  # analyses restrict: being BL stable resident,
  # around for close of BL CHC, resident aged 30+
  restrict <- data.input$resident_0 & data.input$stable_0 &  
    !data.input$ncd_censor_0 & age 
  restrict
}

#*=========
get.pop <- function(data, ncd, outcome){
  
  # who is NCD prevalent at BL or at FUY3
  BL.prev <- FUY3.prev <- rep(0, nrow(data))
  
  if(ncd=='htn'){
    # focus on HTN
    BL.prev[ which(data$htn_0)  ] <- 1
    FUY3.prev[ which(data$htn_3)	 ] <- 1
  } else{
    # hypertensive or DM
    BL.prev[ which(data$dm_0 | data$htn_0) ] <- 1
    FUY3.prev[ which(data$dm_3 | data$htn_3) ] <- 1
  }  
  
  # further restrictions based on the conditioning set for the outcome
  if(outcome=='A'){
    #  does not have any additional restrictions
    restrict <- which(data$U==1)
  } else if(outcome=='B' ){
    # restrict to HIV+ at Y3
    restrict<-  which(data$ncd_hiv_3)
  } else if (outcome=='C'){
    # restrict to NCD-HIV at FUY3
    restrict <- which(FUY3.prev==1 & data$ncd_hiv_3)
  } else if(outcome=='D' ){
    #  prevalent NCD at BL 
    restrict <- which(BL.prev==1)
  } else if (outcome=='E'| outcome=='F') {
    #  prevalent NCD & HIV at BL
    restrict <- which(BL.prev==1 & data$ncd_hiv_0)
  }
  
  data[restrict,]
  
}

preprocess.NCD <- function(data, time, ncd, outcome){
  
  n<- nrow(data)
  
  #	 get measurements for DM
  dm <- get.variables.ncd(data, ncd.string='dm', time=time)
  
  # get measures for HTN	
  htn  <- get.variables.ncd(data, ncd.string='htn', time=time)
  
  # get prevalence, delta, outcome
  if(ncd=='htn'){
    ncd.data <- get.delta.y.single(data=htn, outcome=outcome)
  } else if(ncd=='any'){
    ncd.data <- get.delta.y.anyNCD(data=data, dm=dm, htn=htn, outcome=outcome)
  }
  
  # for the dual-control		
  if(outcome=='C' | outcome=='F'){	
    # note: when used the data are already subsetted on HIV+ status
    TstVL <- supp <- rep(0, n)
    supp.data <- data[, paste('supp', time, sep='_')]
    TstVL[which( !is.na(supp.data ) ) ] <- 1
    supp[ which( supp.data ) ] <- 1 
    # measurement is now both VL and NCD
    delta <- ncd.data$delta*TstVL
    # control of relevant NCD & VL 
    Y <- ncd.data$Y*supp	
  }	else{
    delta <- ncd.data$delta
    Y <- ncd.data$Y
  }
  
  # censoring (no distinction between death/outmigration)
  if(time==0){
    censor <- rep(0, n)
  } else{
    dead <- move <- rep(0, n)
    dead[which(data$dead_3 ) ] <- 1
    move[which(data$outmigrate_3) ] <-1 
    censor <- as.numeric( dead  | move)
  }	
  
  data.frame(censor, delta, prev=ncd.data$prev, Y)
}


# get screening, control, prevalence variables
get.variables.ncd <- function(data, ncd.string, time){
  
  n <- nrow(data)
  
  # set all measurements to 0 
  screen  <-  delta.control <- uncontrol <- control <-  prev <- rep(0, n )
  
  # screening questions
  dx <- data[, paste(ncd.string, 'self_dx', time, sep='_') ]
  txt <- data[, paste(ncd.string, 'self_txt', time, sep='_') ]
  
  screen[ which( !is.na(dx) | !is.na(txt) ) ] <- 1  	
  
  # uncontrol for NCD
  uncontrol.data <- data[, paste(ncd.string, 'uncontrol',time, sep='_')] 

  # set to 0 unless evidence otherwise
  delta.control[ which( !is.na(uncontrol.data)) ] <- 1
  uncontrol[ which( uncontrol.data) ] <- 1
  control[ which( !uncontrol.data) ] <- 1
  
  # prevalent at time t		
  prev[which(data[, paste(ncd.string, time, sep='_') ] ) ] <- 1
  
  #prevalence requires screening & lab measures
  delta.prev <- screen*delta.control
  
  data.frame(delta.control, uncontrol, control, delta.prev, prev)	
}



# get.delta.y.single: indicators for prevalent, measured, control/uncontrol
#   when interestd in a single NCD
get.delta.y.single <- function(data, outcome){
  
  if(outcome=='A' | outcome=='B'){
    # measurement indicator is for prevalence
    delta <- data$delta.prev
    # interested in uncontrolled
    Y<- data$uncontrol
    
  }	else{
    # measurement indicator is for control 
    # (condition on known prevalent at t)
    delta <- data$delta.control
    # interested in controlled
    Y <- data$control	
  }
  data.frame(prev=data$prev, delta, Y)
}


# get.delta.y.any: indicators for prevalent, measured, control/uncontrol
#   when interestd in any NCD
get.delta.y.anyNCD <- function(data, dm, htn, outcome){
  
  # prevalent if have 1+ NCD 
  prev <- as.numeric( htn$prev | dm$prev)
  
  if(outcome=='A' | outcome=='B'){
    # measurement indicator is for prevalence; need screen/labs for both
    delta <- htn$delta.prev*dm$delta.prev			
    # outcome is uncontrolled on either
    Y <-  as.numeric(htn$uncontrol | dm$uncontrol)
  } else{
    
    if(outcome=='D' | outcome=='E' | outcome=='F'){
      # we want control among BL prevalent
      htn.only <- which( data$htn_0 &  (is.na(data$dm_0) | !data$dm_0 ) )
      dm.only <- which( data$dm_0 &  (is.na(data$htn_0) | !data$htn_0 ) )
      both<- which( data$htn_0 & data$dm_0)
    }else{
      # control among FUY3 prevalent	
      htn.only <- which( data$htn_3 &  (is.na(data$dm_3) | !data$dm_3 ) )
      dm.only <- which( data$dm_3 &  (is.na(data$htn_3) | !data$htn_3) )
      both<- which( data$htn_3 & data$dm_3)
    }
    
    # set measurement and actual control to NA
    delta <- Y  <- rep(NA, nrow(data))
    
    # look at measurement among those with the relevant disease
    delta[htn.only] <- htn$delta.control[htn.only]
    delta[dm.only] <- dm$delta.control[dm.only]
    delta[both] <- (htn$delta.control*dm$delta.control)[both]
    
    # look at control among those with the relevant disease
    Y[htn.only] <- htn$control[htn.only]
    Y[dm.only] <- dm$control[dm.only]
    Y[both] <- as.numeric(htn$control*dm$control)[both]
    
  }
  data.frame(prev, delta, Y)
  
}





#*==========
#DO STAGE 1 For NCDS

Stage1.NCD<- function(data.input, settings, outcome, SL.library){
  
  
  # cluster-level variables to retain for stage 2
  E <- c( 'U', 'region_name', 'community_name', 'id', 'pair',  'A', 
          'ncd_control_0', 'dual_control_0',
          'chc_cover_0', 'alcohol_prev_0','overwt_prev_0')
  
  
  # which clusters should be used 
  clusters <- get.clust.exclusion(data.input=data.input, time=settings$time, 
                                  outcome=outcome, 
                                  ncd=settings$ncd)
  nClust <- length(clusters)
  
  stage2 <- data.frame(matrix(NA, nrow=nClust, ncol=(length(E) + 4)) )
  
  if(outcome=='A' | outcome=='B'){
    data.clust <- data.frame(matrix(NA, nrow=nClust, ncol=(length(E) + 18) )) 
  }else{
    data.clust <- data.frame(matrix(NA, nrow=nClust, ncol=(length(E) + 10) )) 
  }
  
  for(j in 1: nClust){
    
    # measured on the whole population of interest		
    these.units <- clusters[j]== data.input$id
    OC <- data.input[these.units, ]
    
    # restrict to the relevant population of interest
    OC <- get.pop(data=OC, ncd=settings$ncd, outcome=outcome)
    
    est <- 	do.NCD.analysis(data=OC, 	
                            outcome=outcome,
                            settings=settings,
                            SL.library= SL.library)
    
    print(j)
    # # create-level dataframe
    data.clust[j,] <- c( OC[1,E], est$data.clust) 
    stage2[j, ] <-  c( OC[1,E], est$stage2)
    
  }
  colnames(data.clust) <- c(E,colnames(est$data.clust))
  colnames(stage2) <- c(E,colnames(est$stage2) )
  
  
  list(data.clust=data.clust, stage2=stage2, W=est$W)
}






# want to exclude clusters that did not measure the relevant NCD at time t
# if BL then looking at prevalence and control 
# if FUY3 then outcomes outcomes D-F restrict to known BL prev.
get.clust.exclusion<- function(data.input, time, outcome, ncd){
  
  clusters <- unique(data.input$id)
  
  # Rely on measured BL NCD status
  if(time==0 | outcome=='D' | outcome=='E' | outcome=='F'){
    
    skip <- rep(F, length(clusters))
    
    for(j in 1: length(clusters)){
      these.units <- clusters[j]== data.input$id
      OC <- data.input[these.units, ]
      
      if(ncd=='any'){
        # if any NCD, require diabetes at BL
        skip[j] <- nrow(OC)==sum(is.na(OC$dm_uncontrol_0 )) 
        
        } else{
          # otherwise only require HTN at BL
          skip[j] <-  nrow(OC)==sum(is.na(OC$htn_uncontrol_0))
        }
      }
    clusters<- clusters[!skip]
  }
  clusters
}  



#*======
do.NCD.analysis<- function(data, outcome,
                           settings, SL.library){

  # first for the control outcomes
  out.control <- do.tmle.control(data=data, outcome=outcome,
                                 settings=settings,
                                 SL.library=SL.library)
  # returns proportion uncontrolled for outcomes A & B
  W <- out.control$W
  
  if(outcome!='A' & outcome!='B'){
    # interested in control among known prevalent
    Ns<- out.control$Ns
    colnames(Ns)<- c('N.pop', 'N.meas', 'N.control', 'N.control.U')
    
    control <- out.control$control.adj$e
    controlU <- out.control$control.unadj$e
    
    stage2<- data.frame(
      Yc= control$pt,
      YcU= controlU$pt,
      nIndv_Yc= Ns$N.pop,
      nIndv_YcU= Ns$N.pop)

    colnames(control)<- paste('control', colnames(control))
    colnames(controlU)<- paste('control.U', colnames(controlU))
    
    data.clust <- data.frame(Ns, control, controlU)
    
  } else{
    # also need to estimate prevalence
    out.prev <- do.tmle.prev(data=data, outcome=outcome, 
                             settings=settings, 
                             SL.library=SL.library)
    # get ratios: proportion uncontrolled, given prevalent
    Bayes.adj <- get.var.bayes(mu1= out.control$control.adj$est$pt, 
                               IC1= out.control$control.adj$IC,
                               mu0= out.prev$prev.adj$est$pt, 
                               IC0= out.prev$prev.adj$IC)$est
    Bayes.unadj <- get.var.bayes(mu1= out.control$control.unadj$est$pt, 
                                 IC1= out.control$control.unadj$IC,
                                 mu0= out.prev$prev.unadj$est$pt,
                                 IC0= out.prev$prev.unadj$IC)$est
    #  this is calculated in terms of uncontrol
    # so flip them
    Bayes.adj <- flip.uncontrol.est(Bayes.adj)
    Bayes.unadj <- flip.uncontrol.est(Bayes.unadj)
    
  
    # outpput
    Ns<- data.frame(out.prev$Ns, out.control$Ns[3:4])
    colnames(Ns)[5:6] <- c('N.uncont', 'N.uncont.U')
    
    stage2<- data.frame(Yc=Bayes.adj$pt, YcU=Bayes.unadj$pt,
                        nIndv_Yc=Ns$N.prev, nIndv_YcU=Ns$N.prev.U)
    # make pretty
    prev<- out.prev$prev.adj$est
    prevU <- out.prev$prev.unadj$est
    colnames(prev)<- paste('prev', colnames(prev))
    colnames(prevU)<- paste('prev.U', colnames(prevU))
    
    colnames(Bayes.adj)<- paste('control', colnames(Bayes.adj))
    colnames(Bayes.unadj)<- paste('control.U', colnames(Bayes.unadj))
    
    data.clust <- data.frame(Ns, prev, prevU, Bayes.adj, Bayes.unadj)
    
    }  
 
   list(stage2=stage2, data.clust=data.clust, W=W)

}

#++++++++++++++++++++++++++++++++++++++++++

# do.tmle.control = 
do.tmle.control <- function(data, outcome, settings, SL.library=NULL){
  
  # get baseline predictors, censoring, measurement, prevalent, control
  baseline.pred <-  get.X(data=data, analysis='NCD', time=settings$time)
  
  # restricted adjustment set to avoid overfitting
  if(outcome=='C' | outcome=='E' | outcome=='F'){
    if(settings$time>0 ){
      baseline.pred <- subset(baseline.pred, 
                              select=c(age.40.49, age.50.59, age.60.plus, 
                                       male, chc.BL))
    }else{
      baseline.pred <- subset(baseline.pred, 
                              select=c(age.40.49, age.50.59, age.60.plus, 
                                       male, mobile))
    }
  }
  # preprocess
  counts <- preprocess.NCD(data=data, time=settings$time, ncd=settings$ncd,
                           outcome=outcome)
  data <- cbind(baseline.pred, counts)
  
  # handle censoring by death or outmigration
  data <- data[!data$censor, ]	
  
  # Number in the population of interest
  N.pop <-  nrow(data)
  
  # setup & run ltmle
  W <- colnames(baseline.pred)
  A<- 'delta' # intervention variable is always delta.
  Y<-'Y'   # Uncontrolled//controlled at time t:
  
  control.adj <- call.ltmle.ncd(data=data, W=W, C=NULL, A=A, Y= Y,
                                SL.library=SL.library)		
  
  N.outcome <- round(control.adj$e$pt*N.pop, 0)
  
  # Secondary analysis
  control.unadj<-  call.ltmle.ncd(data=data, W=NULL, C=NULL, A=A, Y= Y,
                                  SL.library=SL.library)	
  # mean(data[data$delta==1, 'Y'])	
  N.outcome.unadj <- round(control.unadj$e$pt*N.pop, 0)
  
  Ns<- data.frame(N.pop, N.meas=sum(data[,A]),
                  N.out=N.outcome, N.out.unadj=N.outcome.unadj)
  list(Ns=Ns, control.adj=control.adj,
       control.unadj=control.unadj, W=W)
}

call.ltmle.ncd<- function(data, W=NULL, C=NULL,  A,  Y, 
                          SL.library=NULL, 
                          deterministicQ=NULL,
                          observation.weights=NULL, 
                          id=NULL){
  
  
  data.temp<- data[ , c(W, C, A, Y)]
  
  est.temp<- ltmle(data=data.temp, Anodes=A, Cnodes=C, Ynodes=Y, 
                   abar=1,
                   stratify = T, 
                   SL.library=SL.library, 
                   variance.method='ic', 
                   deterministic.Q.function=deterministicQ,
                   observation.weights=observation.weights, 
                   id=id,
                   estimate.time=F)
  
  IC<- est.temp$IC$tmle
  
  est<- data.frame(pt=est.temp$estimate["tmle"], 
                   CI.lo=summary(est.temp)$treatment$CI[1], 
                   CI.hi=summary(est.temp)$treatment$CI[2] 	)
  
  list(est=est,IC=IC)
}


# for outcomes A & B need to estimate prevalence
do.tmle.prev <-  function(data, outcome,
                          settings, SL.library=NULL){
  
  
  # get baseline predictors, censoring, measurement, prevalent, control
  adj <- get.adjustment.prev(data=data, settings=settings)
  counts <- preprocess.NCD(data=data, time=settings$time, ncd=settings$ncd, 
                           outcome=outcome)
  
  data <- cbind(adj$adj, counts)
  
  # handle censoring by death or outmigration
  data <- data[!data$censor, ]	

  # Number in the population of interest
  N.pop <-  nrow(data)
  
  W <- colnames(adj$adj)
  A<- 'delta'
  # outcome as underlying prevalent NCD
  Y<- 'prev'
  
  prev.adj <- call.ltmle.ncd(data=data, W=W, C=NULL, A=A, Y= Y,
                             SL.library=SL.library, deterministicQ=adj$detQ)		
  # Number estimated to be NCD+ = (Estimated prevalence) x (Population size) 
  N.prev <- round(prev.adj$e$pt*N.pop, 0)
  
  #*************** Secondary analysis  ******************************
  prev.unadj<-  call.ltmle.ncd(data=data, W=NULL, C=NULL, A=A, Y= Y,
                               SL.library=SL.library)	
  # sum( data[,A] & data[,Y]) / sum(data[,A])
  N.prev.U <- round(prev.unadj$e$pt*N.pop, 0)
  
  Ns<- data.frame(N.pop, N.meas=sum(data[,A]),
                  N.prev,  N.prev.U=N.prev.U)
  list(Ns=Ns, prev.adj=prev.adj, prev.unadj=prev.unadj, W=W)	
}

get.adjustment.prev <- function(data, settings){
  
  baseline.pred <-  get.X(data=data, analysis='NCD', time=settings$time)
  
  if(settings$time==0){
    # no deterministic knowledge at BL
    adj <- baseline.pred
    detQ <- NULL
  } else { 
    
    # if prevalent at BL, then prevalent at FUY3
    # who is NCD prevalent at BL or at FUY3
    detQ.variable <- rep(0, nrow(data))
   
    if(settings$ncd=='any'){
      detQ.variable[ which(data$dm_0 | data$htn_0) ] <- 1
    } else if (settings$ncd=='htn'){
      detQ.variable[ which(data$htn_0)  ] <- 1
    } 
    
    adj <- data.frame(baseline.pred, detQ.variable)
    
    detQ <- deterministicQ_YES
  }
  list(adj=adj, detQ=detQ)
}
#


flip.uncontrol.est <- function(Bayes){
  pt <- 1 - Bayes$pt
  CI.lo <- 1 - Bayes$CI.hi
  CI.hi <- 1- Bayes$CI.lo
  Bayes.new <- data.frame(pt, CI.lo, CI.hi)
  Bayes.new
}





#*===================================================

Stage2.NCD<- function(data.input, settings, outcome, SL.library){
  
  stage1  <- Stage1.NCD(data.input=data.input, settings=settings, outcome=outcome,
                        SL.library=SL.library)
  
  # primary analysis is for arithmetic risk ratio
  goal <-  'aRR'
  # primary analysis weights indv equally 
  weighting <- 'indv'
  # primary analysis preserves the pairs 
  break.match <- F
  
  clust.adj <-  get.clust.adj(outcome=outcome)
  Yc.unadj <- Stage2(goal=goal, weighting=weighting, data.input=stage1$stage2,
                     outcome= 'YcU', clust.adj=NULL, do.data.adapt=F,  
                     break.match= break.match)
  Yc.adj <- Stage2(goal=goal, weighting=weighting, data.input=stage1$stage2, 
                   outcome='Yc', clust.adj=clust.adj,  do.data.adapt=T,  
                   break.match= break.match)
  list(data.clust=stage1$data.clust, est=rbind(Yc.adj, Yc.unadj), W=stage1$W )
}



# get cluster level adjustment
get.clust.adj<- function(outcome){
  
    if(outcome=='A' | outcome=='B'){
      clust.adj <- c('U', 'chc_cover_0', 'overwt_prev_0')
    } else if(outcome=='C' | outcome=='F'){
      clust.adj <- c('U', 'chc_cover_0', 'dual_control_0')
    } else {
      clust.adj <- c('U', 'chc_cover_0', 'ncd_control_0')
    }
  
  clust.adj
}

#*====== 


get.CI.prev.control<- function(ncd, data.clust, weighting='indv'){
  
  prev <- get.CIs(data.clust=data.clust, weighting=weighting, 
                  nIndv=data.clust$N.pop, Y=data.clust$prev.pt)
  prev.U <- get.CIs(data.clust=data.clust, weighting=weighting, 
                    nIndv=data.clust$N.pop, Y=data.clust$prev.U.pt)
  control <- get.CIs(data.clust=data.clust, weighting=weighting,
                       nIndv=data.clust$N.prev, Y=data.clust$control.pt)
  control.U <- get.CIs(data.clust=data.clust, weighting=weighting, 
                         nIndv=data.clust$N.prev.U, Y=data.clust$control.U.pt)
  
  EST <- data.frame(rbind(prev=unlist(prev), 
                          prev.U=unlist(prev.U),
                          control=unlist(control), 
                          control.U=unlist(control.U) ))
  EST
}

get.CIs <- function(data.clust, weighting, nIndv, Y){
  
  data.temp <- data.frame(id=1:nrow(data.clust), nIndv=nIndv)
  
  # now combine in stage2
  alpha<- get.weights(data.temp, weighting=weighting)$alpha
  
  data.temp<- data.frame(A=rep(1, nrow(data.clust)), Y=Y)
  all <- call.ltmle.ncd(data=data.temp, A='A', Y='Y', 
                         observation.weights=alpha)$est
  colnames(all) <- paste('All', colnames(all), sep='.') 
  
  # among txt
  txt<- data.clust$A==1
  txt <- call.ltmle.ncd(data=data.temp[txt,], A='A', Y='Y', 
                             observation.weights=alpha[txt])$est
  colnames(txt)<- paste('Txt', colnames(txt), sep='.')
  
  con <- data.clust$A==0
  con <- call.ltmle.ncd(data=data.temp[con,], A='A', Y='Y', 
                                 observation.weights=alpha[con])$est
  colnames(con)<- paste('Con', colnames(con), sep='.')
  
  c(all, txt, con)
  
}


#* FILE
get.file.name.ncd <- function(ncd='htn', 
                              time=3, 
                              date=NULL){
  if(is.null(date)){
		date <- 	format(Sys.time(), "%d%b%Y")
	}
 
  file.name <- paste('NCD',  ncd,
                     paste('Yr', time, sep=''),
                     paste('v', date,  sep=''), sep="_")
  
	file.name
}


#===*
  
# MAKE PRETTY OUTPUT

make.pretty.ncd<- function(out){
  out.csv <- data.frame(matrix(NA, nrow=9, ncol=8))
  colnames(out.csv) <- c('', '', '','pt','CI.lo', 'CI.hi', 'pval', 'N.analysis')
  out.csvBL <- out.csv
  
  # Prevalent at FUY3
  out.csv[1:3,] <- get.outcome.output(est=out$EST, this='A: All* with prevalent ncd at Y3; Adj.', 
                                      data=out$outA)
  out.csv[4:6,]<- get.outcome.output(est=out$EST, this='B: HIV+* with prevalent ncd at Y3; Adj.',
                                     data= out$outB)
  out.csv[7:9,]<- get.outcome.output(est=out$EST, this='C: HIV+:Dual control among prevalent ncd at Y3**; Adj.',
                                     data= out$outC)
  
  # Prevalent at BL 
  out.csvBL[1:3,] <- get.outcome.output(est=out$EST, this='D: All with prevalent ncd at BL; Adj.', 
                                        data=out$outD)
  out.csvBL[4:6,]<- get.outcome.output(est=out$EST, this='E: HIV+ and prevalent ncd at BL; Adj.',
                                       data= out$outE)
  out.csvBL[7:9,]<- get.outcome.output(est=out$EST, this='F: HIV+:Dual control among prevalent ncd at BL**; Adj.',
                                       data= out$outF)
  
  out.csv[,3] <- out.csvBL[,3] <- rep(c('Intervention', 'Control', 'RR'),3)
  
  out.csv[,1] <- 'Among prevalent HT at year 3'
  out.csvBL[,1] <- 'Among prevalent HT at baseline'
  out.csv[,2] <- c(rep('All',3), rep('HIV+ at year 3', 3), rep('Dual HTN & HIV RNA control',3) )
  out.csvBL[,2]<- c(rep('All',3), rep('HIV+ at baseline', 3), rep('Dual HTN & HIV RNA control',3) )
  
  OUT <- rbind(out.csv, out.csvBL)
  OUT
}


get.outcome.output <- function(est, this, data){
  print(this)
  out.csv <- data.frame(matrix(NA, nrow=3, ncol=8))
  colnames(out.csv) <- c('', '', '','pt','CI.lo', 'CI.hi', 'pval', 'N.analysis')
  
  txt <- data[data$A==1,]
  con <- data[data$A==0,]
  out.csv[1, c('pt', 'CI.lo', 'CI.hi', 'N.analysis')] <- c(est[this, c('Txt.est', 'Txt.CI.lo', 'Txt.CI.hi')], 
                                                           sum(txt$N.pop) )
  out.csv[2, c('pt', 'CI.lo', 'CI.hi', 'N.analysis')] <- c(est[this, c('Con.est', 'Con.CI.lo', 'Con.CI.hi')], 
                                                           sum(con$N.pop) )
  out.csv[3, c('pt', 'CI.lo', 'CI.hi', 'pval', 'N.analysis')] <- c(est[this, c('Effect.est', 'Effect.CI.lo', 
                                                                               'Effect.CI.hi', 'pval')],
                                                                   sum(data$N.pop))
  out.csv
}
