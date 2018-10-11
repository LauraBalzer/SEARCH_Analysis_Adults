##############
# Master function to  complete the primary & secondary analyses for NCDs
#  detailed in the SEARCH Analysis Plan
#
# Code for additional sensitivity analyses are available by request 
#
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Study Statistician for SEARCH
#---------------------------
#
##############

rm(list=ls())

set.seed(1)

library('SuperLearner')
library('ltmle')

# load functions to do estimation 

source('Stage1_Functions.R')
source('NCD_Functions.R')

source('Stage2_Functions.R')
source('Adapt_Functions.R')


DATE <- 'Final'

# time=3 returns proportion controlled at FUY3
time <- 3


# ncd<- 'htn': with prevalent HT
  # ncd='any': prevalent hypertension or diabetes 
ncd <- 'htn'

# Outputs: Proportion with Control at FUY3
# outcome A: All* with prevalent ncd at FUY3  
# outcome B: HIV+* with prevalent ncd at FUY3
# outcome C: HIV+: Dual control among prevalent ncd at FUY3**  
# outcome D: All with prevalent ncd at BL  
# outcome E: HIV+ and prevalent ncd at BL  
# outcome F: HIV+:Dual control among prevalent ncd at BL**   
#
# all measures adjusted for incomplete measures of ncd control using individual-level TMLE;
# *further adjusted for incomplete measures of ncd prevalence; 
# **further adjusted for incomplete measures of viral suppression.
# Comparison between arms based on community-level TMLE
#########

SL.library=list(
	c('SL.glmnet', 'screen.corP'),
	c('SL.glm', 'screen.corP'),
	c('SL.mean', 'All')
)

## for debugging only, 
## SL.library <- 'glm'

settings<- data.frame(ncd, time)

file.name= get.file.name.ncd(ncd=ncd, time=time, date=DATE)

out<- Run.NCD(settings=settings,  SL.library=SL.library)

save(settings,  out, SL.library, file=paste(file.name, 'Rdata', sep='.'))

write.csv(make.pretty.ncd(out), row.names=F, file=paste(file.name, 'csv', sep='.'))
