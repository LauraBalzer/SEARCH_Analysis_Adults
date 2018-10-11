##############
# R code to complete the analyses for the HIV care cascade
#	as detailed in the SEARCH Analysis Plan
#
# Code for additional sensitivity analyses is available by request 
#
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Study Statistician for SEARCH
###############

rm(list=ls())


set.seed(1)
library('SuperLearner')
library('ltmle')

source('Cascade_Functions.R')
source('Stage1_Functions.R')

source('Stage2_Functions.R')
source('Adapt_Functions.R')


DATE <- 'final'

# specifiy the follow-up year
# options: 0, 1 (intervention only), 2 (intervention only), 3 (outputs comparison)
time <- 3

# primary is fully adjusted 
# (use smaller adj set only for debugging)
adj.full <- T


#========================================================
#########
# overall or region specific
region <- 'All'
#region <- 'EU'
#region <- "SWU"
#region <- 'Kenya'


# within region=='All', prespecified subgroups
# if region!='All', subgroup should = 'All'
subgroup <- "All"
# subgroup <- "Male" 				 
# subgroup <- "Female"
# subgroup <- 'Young'  
# subgroup <- 'Old'

# STAGE1 SL.Library 
SL.library <- list( c('SL.glm', 'screen.corRank10'), 
                    c('SL.glm', 'screen.corP'),
                    c('SL.gam', 'screen.corRank10'), 
                    c('SL.gam', 'screen.corP'),
                    c('SL.mean', 'All'))

# STAGE2 ADJUSTMENT
clust.adj <- c("U","pSupp_0", "pAdol_0")

# SUBGROUP MODIFICATION
if (subgroup=='Young' ){
  # if youth subgroup, then stage1 with glm & limited adjustment set
  SL.library <- 'glm'
  adj.full=F
}
if(region=='Kenya' | region=='EU' | region=='SWU'){
  # if region-specific subgroup, then stage2 with unadjusted
	clust.adj <- 'U'
} 



#########################################################


file.name= get.file.name.cascade(region=region, 
                                 subgroup=subgroup, 
                                 time=time,
                                 date=DATE)

settings<- data.frame(region, subgroup, time, adj.full)	

rm(region, subgroup, time, adj.full)
	
out<- Run.Cascade(settings=settings, SL.library=SL.library, clust.adj= clust.adj)



save(settings, clust.adj, out, SL.library, file=paste(file.name, 'Rdata', sep='.'))



sum(out$data[out$data$A==1, 'N.TstVL'])
sum(out$data[out$data$A==0, 'N.TstVL'])


if(settings$time==3 & settings$region=='All' & settings$subgroup=='All'){
  # reorganize the supplementary table
  
  sec1 <- out$raw2.sec1
  sec0 <- out$raw2.sec0
  rownames(sec1) <- paste(rownames(sec1), 1, sep='.')
  rownames(sec0) <- paste(rownames(sec0), 0, sep='.')
  
  supl.csv <- rbind(
    sec1['pDx.pos.1',],
    sec0['pDx.pos.0',],
    sec1['eART.pDx.1',],
    sec0['eART.pDx.0',],
    sec1['supp.eART.1',],
    sec0['supp.eART.0',],
    sec1['supp.pos.1',],
    sec0['supp.pos.0',],  
    sec1['pos.noTstVL.1',],
    sec0['pos.noTstVL.0',],
    sec1['allpos.noTstVL.1',],
    sec0['allpos.noTstVL.0',]  
  )
  write.csv(supl.csv, file='supl.cascade.final.csv')
}


# # ONCE THE PRIMARY ANALYSIS AND PRESPECIFIED SUBGROUPS HAVE BEEN RUN
# 
# CSV <- rbind(
#   make.pretty.min(time=0),
#   make.pretty.min(time=1),
#   make.pretty.min(time=2),
#   make.pretty.min(time=3),
#   '',
#   make.pretty.min(region='SWU', time=0),
#   make.pretty.min(region='SWU', time=1),
#   make.pretty.min(region='SWU', time=2),
#   make.pretty.min(region='SWU',time=3),
#   '',
#   make.pretty.min(region='EU', time=0),
#   make.pretty.min(region='EU', time=1),
#   make.pretty.min(region='EU', time=2),
#   make.pretty.min(region='EU',time=3),
#   '',
#   make.pretty.min(region='Kenya', time=0),
#   make.pretty.min(region='Kenya', time=1),
#   make.pretty.min(region='Kenya', time=2),
#   make.pretty.min(region='Kenya',time=3),
#   '',
#   make.pretty.min(subgroup='Female',time=3),
#   '',
#   make.pretty.min(subgroup='Male',time=3),
#   '',
#   make.pretty.min(subgroup='Young',time=3),
#   '',
#   make.pretty.min(subgroup='Old',time=3)
# )
# 
# write.csv(CSV, row.names=F, file='Cascade.final.csv')
        