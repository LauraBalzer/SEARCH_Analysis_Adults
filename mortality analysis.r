# By Yea-Hung Chen

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

# sec3pri: load data
load('prepared data.RData')

# sec3pri: define risk period
ii$start<-as.Date('1900-01-01')
is.na(ii$start)<-TRUE
ii$start<-if_else(ii$age_0<15,
                  ii$chc_start_0+(15-ii$age_0)*365.25,
                  ii$start)
is.na(ii$start)<-(ii$start<ii$chc_start_0)
is.na(ii$start)<-(ii$trk_end_3<ii$start)
ii$start<-pmax(ii$start,ii$chc_start_0,na.rm=TRUE)
ii$end<-ii$trk_end_3

# sec3pri: impute missing dates
ii<-impute_dates(c('di.end','do.end','om.end'))

# sec3pri: define survival variables
ovar<-c('di')
cvar<-c('do','om')
ii<-define_survival_variables(ii)

# sec3pri: define population
ii<-subset(ii,!data_flag)
ii<-subset(ii,!(dead_0|move_0))
ii<-subset(ii,!is.na(dead_3))
ii<-subset(ii,resident_0)
ii<-subset(ii,stable_0)

# sec3pri: define adjustment variables
ii$hiv<-ifelse(ii$hiv_0,1,0)
ii$wealth<-ifelse(ii$wealth_0=='0',1,0)

# sec3pri: initiate community-level data frame
cc<-c('com','community_name','pair','intervention')
cc<-unique(ii[,cc])
cc<-cc[order(cc$com),]
row.names(cc)<-NULL  

# sec3pri: outcome variables
cc<-oo('sec3pri',method='rate')

# sec3pri: add adjustment variables to community-level data frame
aa('hiv')
aa('wealth')

# sec3pri: pre-process data
cc<-preprocess(cc,YHC=TRUE)

# sec3pri: stage-ii analysis
s2('sec3pri',CLUST.ADJ=c('U','hiv','wealth'),SURVIVAL=FALSE)

# sec1pri: load data
load('prepared data.RData')

# sec1pri: define risk period
ii$start<-ii$chc_start_0
ii$end<-ii$trk_end_3

# sec1pri: impute missing dates
ii<-impute_dates(c('di.end','do.end','om.end'))

# sec1pri: define survival variables
ovar<-c('di')
cvar<-c('do','om')
ii<-define_survival_variables(ii)

# sec1pri: define population
ii<-subset(ii,!data_flag)
ii<-subset(ii,!(dead_0|move_0))
ii<-subset(ii,!is.na(dead_3))
ii<-subset(ii,adult_0)
ii<-subset(ii,resident_0)
ii<-subset(ii,hiv_0)
ii<-subset(ii,stable_0)

# sec1pri: define adjustment variables
ii$cd4.350<-with(ii,ifelse(cd4_0<=350,1,0))
ii$cd4.50<-with(ii,ifelse(cd4_0<=50,1,0))
ii$vl.500<-ifelse(ii$supp_0,1,0)

# sec1pri: initiate community-level data frame
cc<-c('com','community_name','pair','intervention')
cc<-unique(ii[,cc])
cc<-cc[order(cc$com),]
row.names(cc)<-NULL  

# sec1pri: outcome variables 
cc<-oo('sec1pri',method='survival')

# sec1pri: add adjustment variables to community-level data frame
aa('cd4.350')
aa('cd4.50')
aa('vl.500')

# sec1pri: pre-process data
cc<-preprocess(cc,YHC=TRUE)

# sec1pri: stage-ii analysis
s2('sec1pri',CLUST.ADJ=c('U','cd4.350','cd4.50','vl.500'),SURVIVAL=TRUE)

# output estimates
tt<-rbind(sec3pri,sec1pri)
write.csv(tt,'mortality intervention effects.csv',row.names=FALSE)