# Stage1Functions_.R 
# Helper functions for Stage 1: estimating the community-level outcomes & getting covariate data
#
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Lead Statistician for SEARCH



get.subgroup <- function(data, subgroup, time=0){
	
  subgroup <- as.character(subgroup)
	this.subgroup <- rep(F, nrow(data))

	youth <- get.age.grp(data, time)
	if(is.null(subgroup) ){
		subgroup <- 'All'
	}
	
	if(subgroup=='All'){
		# if no subgroups of interest
		this.subgroup[1:nrow(data) ] <- T
	
	} else if (subgroup=='EU'){
		this.subgroup[ which(data$region_name=='Eastern Uganda') ] <- T  
	} else if (subgroup=='SWU'){
		this.subgroup[ which(data$region_name=='Western Uganda') ] <- T  
	} else if (subgroup=='Kenya'){
		this.subgroup[ which(data$region_name=='Kenya') ] <- T  
	# SEX
	} else if(subgroup=='Male'){
		this.subgroup[ which(data$sex_0) ] <- T
	} else if(subgroup=='Female'){
		this.subgroup[ which(!data$sex_0) ] <- T
	# AGE
	} else if(subgroup=='Young'){
		this.subgroup[ youth ] <- T
	} else if(subgroup=='Old'){
		this.subgroup[ !youth] <- T
	# MOBILITY AND VMC	
	} else if(subgroup=='NonMobile'){
		this.subgroup[ which(data$moAway_0 <  1) ] <- T
	} else if(subgroup=='UncircMen'){
		this.subgroup[ which(data$non_circum_0) ] <- T
	} 
	print(c(time, subgroup))
	this.subgroup	
}

get.age.grp<- function(data, time=0){
  
  if(time==0){
    youth <- data$age_0 < 25
  } else if(time==1){
    youth <- data$age_0 < 24
  } else if(time==2){
    youth <- data$age_0 < 23
  } else{
    youth <- data$age_0 < 22
  }
  youth
}


# get relevant covariates for predicting Delta & Censoring
get.X <- function(data, analysis='HIV', time=3, adj.full=T){
			
	n <- nrow(data)
	
	# age # reference age group <20
	age.20.29 <- age.30.39 <- age.40.49 <- age.50.59 <- age.60.plus <- rep(0, n)
	age.20.29[ which(data$age_0>19 & data$age_0<30) ] <- 1
	age.30.39[ which(data$age_0>29 & data$age_0<40) ] <- 1
	age.40.49[ which(data$age_0>39 & data$age_0<50) ] <- 1
	age.50.59 [which(data$age_0>49 & data$age_0<60) ] <- 1
	age.60.plus[which(data$age_0>59)  ] <- 1
	age.matrix <- data.frame(cbind(age.20.29, age.30.39, age.40.49, age.50.59, age.60.plus))
	
	# reference is missing
	single <- married <- widowed <- 	divorced.separated <- rep(0, n)
	single[ which(data$marital_0==1)] <- 1
	married[ which(data$marital_0==2) ] <-1
	widowed[ which(data$marital_0 ==3)] <-1
	divorced.separated[ which(data$marital_0==4 | data$marital_0==5)] <-1
	marital <- data.frame(single, married, widowed, divorced.separated)
		
	# education: reference is less than primary or missing
	primary <- as.numeric(data$edu_primary_0)
	secondary.plus <- as.numeric(data$edu_secondary_plus_0)
	education <- data.frame(primary, secondary.plus)

	# occupation: reference NA
	formal.hi <- as.numeric(data$formal_hi_occup_0)
	informal.hi <- as.numeric(data$informal_hi_occup_0)
	informal.lo <- as.numeric(data$informal_low_occup_0)
	jobless <- as.numeric(data$jobless_0)
	student <- as.numeric(data$student_0)
	fisherman <- as.numeric(data$fisherman_0)
	occupation<- data.frame(formal.hi, informal.hi, informal.lo, jobless, student, fisherman)	
	
	# alcohol use: ref is NA
	alcohol.yes <- alcohol.no <- rep(0, n)
	alcohol.yes[which(data$alcohol_0) ] <- 1
	alcohol.no[which(!data$alcohol_0) ] <- 1

	# reference wealth is NA missing
	wealth0 <- wealth1<- wealth2 <- wealth3 <- wealth4 <-  rep(0, n) 
	wealth0[ which(data$wealth_0==0)] <- 1
	wealth1[ which(data$wealth_0==1)] <- 1
	wealth2[ which(data$wealth_0==2)] <- 1
	wealth3[ which(data$wealth_0==3)] <- 1
	wealth4[ which(data$wealth_0==4)] <- 1
	wealth <- data.frame(cbind(wealth0, wealth1, wealth2, wealth3, wealth4))

	#mobility indicators
	mobile <- as.numeric(data$mobile_0)
	
	# shifted main residence
	shift.no <- shift.yes <- rep(0,n)
	shift.no[which(!data$shifted_0)] <-1
	shift.yes[which(data$shifted_0)] <-1
	
	# nights home
	nights <- as.numeric(as.character(data$nightsHome_0))
	nights0 <- nights1.2 <- nights3.4 <- nights5 <- rep(0,n)
	nights0[which(nights==0)] <-1
	nights1.2[which(nights==1 | nights==2)] <-1
	nights3.4[which(nights==3 | nights==4)] <-1
	nights5[which(nights==5)] <- 1
	mobility <- data.frame(mobile, shift.no, shift.yes, nights0, nights1.2, nights3.4, nights5)
	
	# health-seeking
	chc.BL <- as.numeric(data$chc_0)
	self.hivtest.yes <- self.hivtest.no <- rep(0,n)
	self.hivtest.yes[which(data$self_hivtest_0)]<-1
	self.hivtest.no[which(!data$self_hivtest_0)] <-1
	health<- data.frame(chc.BL, self.hivtest.yes, self.hivtest.no)

	male <- rep(0,n)
	male[which(data$sex_0)] <- 1

	X <- cbind(
		age.matrix, marital,
		education, occupation,
		alcohol.yes, alcohol.no, wealth, mobility, male)
		
	if(analysis=='HIV'){
		X<- cbind(X, health)
	
	} else if(analysis=='NCD'){
	  # reference is underweight or NA
	  #set NA if <15 or >40
	  bmi <- data$bmi_0
	  bmi[which(bmi<15)] <- NA
	  bmi[which(bmi>40)] <- NA
	  bmi.norm <- bmi.over <- bmi.obese <- rep(0,n)
	  bmi.norm[ which(bmi >=18 & bmi <25) ] <-1
	  bmi.over[ which(bmi >=25 & bmi <30) ] <-1
	  bmi.obese[ which(bmi >= 30) ] <-1
	  X<- cbind(X, bmi.norm, bmi.over, bmi.obese )
	  X<- subset(X, select=- c( age.20.29,age.30.39, alcohol.no) )
	  
	  if(time>0){
	    # adjust for baseline CHC attendance
	    X<- cbind(X, chc.BL)
	  }
	} else if(analysis=='Cascade' & !adj.full){
	    X <- data.frame(cbind(mobile, male))
	 
	}
	X
	
}




# get.var - function to get inference via the delta method
# 		assumes inputed estimators are asymptotically linear
#		i.e. written in first order as an empircal mean of an influence curve (IC)
#	input:  point estimates (mu1, mu0), corresponding influence curves (IC1, IC0)
#		significance level
#	output: point estimate, var, wald-type CI 

get.var.bayes <- function(mu1, mu0=NULL, IC1, IC0=NULL, alpha=0.05){
	
	mu1<- unlist(mu1)
		
	if(is.null(mu0)){ 
		# if single TMLE 
		psi<- mu1
		IC<- IC1
		log= F
	
	} else { 
		# if ratio of TMLEs (i.e. target = psi/psi0)
		mu0<- unlist(mu0)
		# get inference via the delta method on log scale
		psi<- log(mu1/mu0)
		IC <- 1/mu1*(IC1) - 1/mu0*IC0
		log=T
	}
	
	# variance of asy lin est is var(IC)/n
	var<- var(IC)/length(IC)
	# testing and CI	
	cutoff <- qnorm(alpha/2, lower.tail=F)
	se<- sqrt(var)
	CI.lo <- psi - cutoff*se
	CI.hi <- psi + cutoff*se

	if(log){
		est<- data.frame(pt=exp(psi), CI.lo=exp(CI.lo), CI.hi=exp(CI.hi) ) 
	}else{
		est<- data.frame(pt=psi, CI.lo=CI.lo, CI.hi=CI.hi) 
	}

	list(est=est, IC=IC)
}




#===================================================#===================================================
# SCREENING ALGORITHMS FOR SUPERLEARNER
# See SuperLearner help file for more info: ?SuperLearner
#===================================================#===================================================
screen.corRank10 <- function(Y, X, family, ...) screen.corRank(Y, X, family, rank = 10, ...)
screen.corRank20 <- function(Y, X, family, ...) screen.corRank(Y, X, family, rank = 20, ...)
screen.corRank5 <- function(Y, X, family, ...) screen.corRank(Y, X, family, rank = 5, ...)
screen.corRank3 <- function(Y, X, family, ...) screen.corRank(Y, X, family, rank = 3, ...)

screen.corP3<- function(Y, X, family, ...) screen.corP(Y, X, family, minscreen = 3, ...)

#===================================================#===================================================
# FUNCTIONS TO ENCODE OUR DETERMINISTIC KNOWLEDGE
# See ltmle help file for more info: ?ltmle
# Also see the Analysis Plan
#===================================================#===================================================

# deterministicQ_YES
# if detQ.variable==1, then outcome==1 with probability 1
deterministicQ_YES<- function(data, current.node, nodes, called.from.estimate.g) {
  L2.index <- which(names(data) == "detQ.variable")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data[,L2.index]==1
  return(list(is.deterministic=is.deterministic, Q.value=1))
}

# deterministicQ_NO
# if detQ.variable==0,  then outcome==0 with probability 1
deterministicQ_NO<- function(data, current.node, nodes, called.from.estimate.g) {
  L2.index <- which(names(data) == "detQ.variable")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data[,L2.index]== 0
  return(list(is.deterministic=is.deterministic, Q.value=0))
}

# deterministicQ_combo 
# cannot be suppressed if dead, outmigrated or not on ART
# cannot have Z*=1 if combo= (D=1 OR M=1 OR eART=0)
deterministicQ_combo<- function(data, current.node, nodes, called.from.estimate.g) {
  L2.index <- which(names(data) == "combo")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data[,L2.index]==1
  return(list(is.deterministic=is.deterministic, Q.value=0))
}



