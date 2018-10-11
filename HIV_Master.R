##############
# R code to complete the analyses for the HIV primary & secondary outcomes 
#	detailed in the SEARCH Analysis Plan
#
# Code for additional sensitivity analyses is available by request 
#
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Study Statistician for SEARCH
##############

rm(list=ls())

set.seed(1)
library('SuperLearner')

library('ltmle')


source('HIV_Functions.R')
source('Stage1_Functions.R')

source('Stage2_Functions.R')
source('Adapt_Functions.R')

date <- 'final'

# primary outcome specifications used throughout.
# see analysis plan for secondary/sensitivity parameterizations


#========================================================
#########
# primary analysis includes all individuals
subgroup <- 'All';
#subgroup <- 'EU'
#subgroup <- "SWU"
#subgroup <- 'Kenya'
#subgroup <- "Male" 				 
#subgroup <- "Female"
#subgroup <- 'Young'  
#subgroup <- 'Old'
#subgroup <- 'NonMobile'
#subgroup <- 'UncircMen'



# STAGE1 SL.Library for C/Delta
SL.library=list(
	c('SL.glmnet', 'screen.corP'),
	c('SL.glm', 'screen.corP'),
	c('SL.mean', 'All')
)

#########################################################


out<- Run.HIV(subgroup=subgroup, SL.library=SL.library)
file.name= get.file.name.hiv(subgroup=subgroup, SL.library=SL.library, date=date)

save(out, SL.library, file=paste(file.name, 'Rdata', sep='.'))


#--------

# Once have run the primary & prespecified subgrups,
# output the results for the primary analytic approach

# date <- 'final'

# CSV <- rbind(
#   make.pretty.hiv.min(subgroup='All', date),
#   make.pretty.hiv.min(subgroup='Kenya', date),
#   make.pretty.hiv.min(subgroup='SWU', date),
#   make.pretty.hiv.min(subgroup='EU',date),
#   make.pretty.hiv.min(subgroup='Male', date),
#   make.pretty.hiv.min(subgroup='Female', date),
#   make.pretty.hiv.min(subgroup='Young', date),
#   make.pretty.hiv.min(subgroup='Old', date),
#   make.pretty.hiv.min(subgroup='NonMobile',  date),
#   make.pretty.hiv.min(subgroup='UncircMen',  date))
# write.csv(CSV, row.names=F, file=paste('Primary', date, 'csv', sep='.'))
