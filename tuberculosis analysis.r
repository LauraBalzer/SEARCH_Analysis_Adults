# By Yea-Hung Chen

# clear work space
rm(list=ls())

# set random seed
set.seed(2018) 

# define time point of interest for survival analysis
TIME<-3

# load libraries
library(survival)

# load functions
source('impute dates.r')
source('define survival variables.r')
source('community-level estimates.r')
source('Preprocess_Functions.R')
source('Stage2_Functions.R')
source('Adapt_Functions.R')
source('s2.r')

# pri: load data
load('prepared data.RData')

# pri: define risk period
ii$start<-ii$tb_risk_start
ii$end<-ii$tb_risk_end

# pri: impute missing dates
ii<-impute_dates(c('ht.end','di.end','do.end','om.end'))

# pri: define survival variables
ovar<-c('ht','di')
cvar<-c('do','om')
ii<-define_survival_variables(ii)

# pri: define population
ii<-subset(ii,!data_flag)
ii<-subset(ii,!tb_censor_0)
ii<-subset(ii,!(dead_0|move_0))
ii<-subset(ii,resident_0)
ii<-subset(ii,adult_0)
ii<-subset(ii,!tb_0)
ii<-subset(ii,hiv_0|is.na(hiv_0))
ii<-subset(ii,stable_0)

# pri: initiate community-level data frame
cc<-c('com','community_name','pair','intervention',
      'tb_incidence_0','hiv_prev_0')
cc<-unique(ii[,cc])
cc<-cc[order(cc$com),]
row.names(cc)<-NULL  

# pri: outcome variables
cc<-oo('pri',method='survival')
cc<-oo('pri_cd4_1',method='survival',dd=subset(ii,hiv_0&(cd4_0<=500)))

# pri: pre-process data
cc<-preprocess(cc,YHC=TRUE)

# pri: stage-ii analysis
s2('pri',CLUST.ADJ=c('U','tb_incidence_0','hiv_prev_0'),SURVIVAL=TRUE)
s2('pri_cd4_1',CLUST.ADJ=c('U','tb_incidence_0','hiv_prev_0'),SURVIVAL=TRUE)

# annual: load data
load('prepared data.RData')

# annual: define risk periods
source('risk periods for annual rates.r')

# annual: impute missing dates
ii<-impute_dates(c('tb.end','da.end','om.end'))

# annual: define survival variables
ovar<-c('tb')
cvar<-c('da','om')
ii<-define_survival_variables(ii)

# annual: define population
ii<-subset(ii,!data_flag)
ii<-subset(ii,!tb_censor_0)
ii<-subset(ii,!(dead_0|move_0))
ii<-subset(ii,adult_0)
ii<-subset(ii,resident_0)
ii<-subset(ii,!tb_0)

# annual: initiate community-level data frame
cc<-c('com','community_name','pair','intervention',
      'tb_incidence_0','hiv_prev_0')
cc<-unique(ii[,cc])
cc<-cc[order(cc$com),]
row.names(cc)<-NULL  

# annual: outcome variables
oo('neg',method='annual',dd=subset(ii,!hiv_0))
oo('pos',method='annual',dd=subset(ii,hiv_0))

# annual: pre-process data
cc<-preprocess(cc,YHC=TRUE)

# annual: stage-ii analysis
OUTCOME<-paste(c('neg','pos'),rep(1:3,each=2),sep='')
invisible(lapply(OUTCOME,s2,CLUST.ADJ=c('U','tb_incidence_0','hiv_prev_0'),
                 SURVIVAL=FALSE))

# output estimates
tt<-rbind(pri,pri_cd4_1,neg1,neg2,neg3,pos1,pos2,pos3)
write.csv(tt,'tuberculosis results.csv',row.names=FALSE)