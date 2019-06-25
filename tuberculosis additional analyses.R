# By Yea-Hung Chen

# Restrict to HIV-positive individuals
# Includes analysis for both survival and rates (adjusted/unadjusted)


#------------------------------------------------------
# OUTCOME: HIV-RELATED TUBERCULOSIS OR DEATH FROM ILLNESS 
#------------------------------------------------------

# clear work space
rm(list=ls())

# set random seed
set.seed(2018) 

# define time point of interest for survival analysis
TIME<-3

# load libraries
library(survival)
if_else<-dplyr::if_else

# load functions
source('impute dates.r')
source('define survival variables.r')
source('community-level estimates.r')
source('Preprocess_Functions.R')
source('Stage2_Functions.R')
source('Adapt_Functions.R')
source('s2.r')

# load data
load('prepared data.RData')

# define risk period
ii$start<-ii$tb_risk_start
ii$end<-ii$tb_risk_end

# impute missing dates
ii<-impute_dates(c('ht.end','di.end','do.end','om.end'))

# define survival variables
ovar<-c('ht','di')
cvar<-c('do','om')
ii<-define_survival_variables(ii)

# define population
ii<-subset(ii,!data_flag)
ii<-subset(ii,!tb_censor_0)
ii<-subset(ii,!(dead_0|move_0))
ii<-subset(ii,adult_0)
ii<-subset(ii,resident_0)
ii<-subset(ii,!tb_0)
ii<-subset(ii,hiv_0|is.na(hiv_0))
ii<-subset(ii,stable_0)
ii<-subset(ii,hiv_0) # restricted to known HIV+

# initiate community-level data frame
cc<-c('com','community_name','pair','intervention',
      'tb_incidence_0','hiv_prev_0')
cc<-unique(ii[,cc])
cc<-cc[order(cc$com),]
row.names(cc)<-NULL  

# backup initial community-level data frame
cc.backup<-cc

# risks: outcome variables
cc<-oo('pri',method='survival')
cc<-oo('pri_cd4_1',method='survival',dd=subset(ii,hiv_0&(cd4_0<=500)))

# risks: pre-process data
cc<-preprocess(cc,YHC=TRUE)

# risks: stage-ii analysis
s2('pri',CLUST.ADJ='U',SURVIVAL=TRUE)
s2('pri_cd4_1',CLUST.ADJ='U',SURVIVAL=TRUE)
risk_unadj<-pri
risk_cd4_unadj<-pri_cd4_1
s2('pri',CLUST.ADJ=c('U','tb_incidence_0','hiv_prev_0'),SURVIVAL=TRUE)
s2('pri_cd4_1',CLUST.ADJ=c('U','tb_incidence_0','hiv_prev_0'),SURVIVAL=TRUE)
risk_adj<-pri
risk_cd4_adj<-pri_cd4_1

# rates: outcome variables
cc<-cc.backup
cc<-oo('pri',method='rate')
cc<-oo('pri_cd4_1',method='rate',dd=subset(ii,hiv_0&(cd4_0<=500)))

# rates: pre-process data
cc<-preprocess(cc,YHC=TRUE)

# rates: stage-ii analysis
s2('pri',CLUST.ADJ='U',SURVIVAL=FALSE)
s2('pri_cd4_1',CLUST.ADJ='U',SURVIVAL=FALSE)
rate_unadj<-pri
rate_cd4_unadj<-pri_cd4_1

# collect estimates
tt<-rbind(risk_unadj,risk_adj,rate_unadj,
          risk_cd4_unadj,risk_cd4_adj,rate_cd4_unadj)
tt$measure<-c('unadjusted risk','adjusted risk','unadjusted rate')
tt$population<-c(rep('HIV positive',3),rep('Low CD4',3))
tt<-tt[,c('population','measure',
          'Txt.est','Txt.CI.lo','Txt.CI.hi',
          'Con.est','Con.CI.lo','Con.CI.hi',
          'Effect.est','Effect.CI.lo','Effect.CI.hi')]

# output estimates
write.csv(tt,'primary tuberculosis analysis, restricted.csv',
          row.names=FALSE)

#------------------------------------------------------
# OUTCOME AS HIV-RELATED TUBERCULOSIS
#------------------------------------------------------

# clear work space
rm(list=ls())

# set random seed
set.seed(2018) 

# define time point of interest for survival analysis
TIME<-3

# load libraries
library(survival)
if_else<-dplyr::if_else

# load functions
source('impute dates.r')
source('define survival variables.r')
source('community-level estimates.r')
source('Preprocess_Functions.R')
source('Stage2_Functions.R')
source('Adapt_Functions.R')
source('s2.r')

# load data
load('prepared data.RData')

# define risk period
ii$start<-ii$tb_risk_start
ii$end<-ii$tb_risk_end

# impute missing dates
ii<-impute_dates(c('ht.end','da.end','om.end'))

# define survival variables
ovar<-c('ht')
cvar<-c('da','om')
ii<-define_survival_variables(ii)

# define population
ii<-subset(ii,!data_flag)
ii<-subset(ii,!tb_censor_0)
ii<-subset(ii,!(dead_0|move_0))
ii<-subset(ii,adult_0)
ii<-subset(ii,resident_0)
ii<-subset(ii,!tb_0)
ii<-subset(ii,stable_0)
ii<-subset(ii,hiv_0) # restricted to known HIV+

# initiate community-level data frame
cc<-c('com','community_name','pair','intervention',
      'tb_incidence_0','hiv_prev_0')
cc<-unique(ii[,cc])
cc<-cc[order(cc$com),]
row.names(cc)<-NULL  

# backup initial community-level data frame
cc.backup<-cc

# risks: outcome variables
cc<-oo('sec')
cc<-oo('sec_cd4_1',dd=subset(ii,hiv_0&(cd4_0<=500)))

# risks: pre-process data
cc<-preprocess(cc,YHC=TRUE)

# risks: stage-ii analysis
s2('sec',CLUST.ADJ='U',SURVIVAL=TRUE)
s2('sec_cd4_1',CLUST.ADJ='U',SURVIVAL=TRUE)
risk_unadj<-sec
risk_cd4_unadj<-sec_cd4_1
s2('sec',CLUST.ADJ=c('U','tb_incidence_0','hiv_prev_0'),SURVIVAL=TRUE)
s2('sec_cd4_1',CLUST.ADJ=c('U','tb_incidence_0','hiv_prev_0'),SURVIVAL=TRUE)
risk_adj<-sec
risk_cd4_adj<-sec_cd4_1

# rates: outcome variables
cc<-cc.backup
cc<-oo('sec',method='rate')
cc<-oo('sec_cd4_1',method='rate',dd=subset(ii,hiv_0&(cd4_0<=500)))

# rates: pre-process data
cc<-preprocess(cc,YHC=TRUE)

# rates: stage-ii analysis
s2('sec',CLUST.ADJ='U',SURVIVAL=FALSE)
s2('sec_cd4_1',CLUST.ADJ='U',SURVIVAL=FALSE)
rate_unadj<-sec
rate_cd4_unadj<-sec_cd4_1

# collect estimates
tt<-rbind(risk_unadj,risk_adj,rate_unadj,
          risk_cd4_unadj,risk_cd4_adj,rate_cd4_unadj)
tt$measure<-c('unadjusted risk','adjusted risk','unadjusted rate')
tt$population<-c(rep('HIV positive',3),rep('Low CD4',3))
tt<-tt[,c('population','measure',
          'Txt.est','Txt.CI.lo','Txt.CI.hi',
          'Con.est','Con.CI.lo','Con.CI.hi',
          'Effect.est','Effect.CI.lo','Effect.CI.hi')]

# output estimates
write.csv(tt,'secondary tuberculosis analysis, restricted.csv',
          row.names=FALSE)