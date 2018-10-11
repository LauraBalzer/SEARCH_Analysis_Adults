# By Yea-Hung Chen

# clear work space
rm(list=ls())

# load data
load('outputs-postbl-unblinded.RData')
ii<-data.frame(outputs)

# community
ii$com<-ii$community_number

# active tuberculosis
ii$tb<-ii$tb_3
ii$tb.end<-ii$tb_date
is.na(ii$tb.end)<-!ii$tb

# hiv-associated tuberculosis
ii$ht<-ii$tb_hivasc_3
ii$ht.end<-ii$tb_date
is.na(ii$ht.end)<-!ii$ht

# death
ii$da<-!is.na(ii$dead_3)&ii$dead_3
ii$da.end<-ii$tb_death_date
is.na(ii$da.end)<-!ii$da

# death from illness
ii$di<-with(ii,!is.na(dead_3)&dead_3&!is.na(death_cause_3)&(death_cause_3==1))
ii$di.end<-ii$tb_death_date
is.na(ii$di.end)<-!ii$di

# death from other causes
ii$do<-with(ii,!is.na(dead_3)&dead_3&!is.na(death_cause_3)&(death_cause_3!=1))
ii$do.end<-ii$tb_death_date
is.na(ii$do.end)<-!ii$do

# out-migration
ii$om<-!is.na(ii$outmigrate_3)&ii$outmigrate_3
ii$om.end<-ii$tb_outmigration_date
is.na(ii$om.end)<-!ii$om

# save data
save.image('prepared data.RData')