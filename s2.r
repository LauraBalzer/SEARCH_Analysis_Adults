
# By Yea-Hung Chen

# function for calling Stage2()
s2<-function(OUTCOME,CLUST.ADJ,SURVIVAL){
  ee<-Stage2(outcome=OUTCOME,clust.adj=CLUST.ADJ,
             data.input=cc,survival=SURVIVAL,weighting='indv',
             do.data.adapt=TRUE)
  ee$outcome<-OUTCOME
  ee$adj<-paste(CLUST.ADJ,collapse='+')
  ee$survival<-SURVIVAL
  if(SURVIVAL){ee$time<-TIME} else {ee$time<-'NA'}
  assign(OUTCOME,ee,envir=.GlobalEnv)
}