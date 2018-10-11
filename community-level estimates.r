# By Yea-Hung Chen

# function for adjustment variables
aa<-function(adj,adj.name=NULL,dd=ii){
  # estimates
  aa<-tapply(dd[,adj],dd$com,function(vv){prop.table(table(vv))['1']})
  aa[is.na(aa)]<-0 # applies if no values are equal to 1
  aa<-data.frame(com=as.numeric(names(aa)),adj=as.numeric(aa))
  if(is.null(adj.name)){adj.name<-adj}
  names(aa)<-gsub('adj',adj.name,names(aa))
  # merge estimates into cc
  cc<-merge(cc,aa,by='com',all.x=TRUE,all.y=FALSE)
  # output estimates
  assign('cc',cc,envir=.GlobalEnv)
}

# function for outcome variables
oo<-function(outcome,method='survival',dd=ii){
  nIndv<-paste('nIndv',outcome,sep='_')
  if(method=='survival'){
    cc[,outcome]<-sapply(cc$com,ss,dd=dd)
    cc[,nIndv]<-sapply(cc$com,function(cc){sum(dd$com==cc)})
  } 
  if(method=='rate'){
    cc[,outcome]<-sapply(cc$com,rr,dd=dd)
    cc[,nIndv]<-tapply(dd$rend,dd$com,sum)
  }
  if(method=='annual'){
    cc<-merge(cc,do.call(rbind,lapply(cc$com,yy,cc.var='com',
                                      dd=dd,nn=outcome)),by='com')
  }
  assign('cc',cc,envir=.GlobalEnv)
}

# function for survival
ss<-function(cc,dd,time=TIME){
  dd<-subset(dd,com==cc)
  ss<-with(dd,Surv(time=rend,event=o))
  ee<-survfit(ss~1)
  summary(ee,min(time,max(dd$rend)))$surv
}  

# function for rates
rr<-function(cc,dd){
  dd<-subset(dd,com==cc)
  sum(dd$oend<=dd$rend,na.rm=TRUE)/sum(dd$rend,na.rm=TRUE)
}

# function for annual rates
yy<-function(cc,cc.var='com',dd,nn){
  # subset data
  dd<-dd[dd[,cc.var]==cc,]
  # convert survival end dates back to dates
  dd$oend<-dd$start+dd$oend*365.25
  dd$rend<-dd$start+dd$rend*365.25
  # function for incidence
  incidence<-function(start,end){
    # redefine annual risk period
    is.na(dd$rend)<-!is.na(start)&(dd$rend<start)
    end<-pmin(dd$rend,end,na.rm=TRUE) 
    # outcome
    oo<-!is.na(dd$oend)&!is.na(start)&!is.na(end)&
      (start<=dd$oend)&(dd$oend<=end)
    # person-time
    tt<-as.numeric(end-start)/365.25
    tt<-ifelse(is.na(tt),0,tt)
    # marker for ending risk within annual risk period
    ee<-(start<=dd$rend)&(dd$rend<=end)
    # output results
    list(oo,tt,ee)
  }
  # year-1 incidence rate
  rr1<-incidence(dd$y1start,dd$y1end)
  oo1<-rr1[[1]]
  tt1<-rr1[[2]]
  ee1<-rr1[[3]]
  is.na(dd$y2start)<-ee1
  is.na(dd$y3start)<-ee1
  is.na(dd$y2end)<-ee1
  is.na(dd$y3end)<-ee1
  # year-2 incidence rate
  rr2<-incidence(dd$y2start,dd$y2end)
  oo2<-rr2[[1]]
  tt2<-rr2[[2]]
  ee2<-rr2[[3]]
  is.na(dd$y3start)<-ee2
  is.na(dd$y3end)<-ee2
  # year-3 incidence rate
  rr3<-incidence(dd$y3start,dd$y3end)
  oo3<-rr3[[1]]
  tt3<-rr3[[2]]
  # entire risk period
  dd$start.min<-pmin(dd$y1start,dd$y2start,dd$y3start,na.rm=TRUE)
  dd$end.max<-pmax(dd$y1end,dd$y2end,dd$y3end,na.rm=TRUE)
  rra<-incidence(dd$start.min,dd$end.max)
  ooa<-rra[[1]]
  tta<-rra[[2]]
  # incidence
  rr1<-sum(oo1)/sum(tt1)
  rr2<-sum(oo2)/sum(tt2)
  rr3<-sum(oo3)/sum(tt3)
  # output incidence rates 
  sumoo<-sum(oo1)+sum(oo2)+sum(oo3)
  sumtt<-sum(tt1)+sum(tt2)+sum(tt3)
  if((sum(ooa)==sumoo)&(floor(sum(tta))==floor(sumtt))){
    cc<-data.frame(var=cc)
    names(cc)<-gsub('var',cc.var,names(cc))
    cc[,paste(nn,1,sep='')]<-sum(rr1)
    cc[,paste(nn,2,sep='')]<-sum(rr2)
    cc[,paste(nn,3,sep='')]<-sum(rr3)
    cc[,paste('nIndv_',nn,1,sep='')]<-sum(tt1)
    cc[,paste('nIndv_',nn,2,sep='')]<-sum(tt2)
    cc[,paste('nIndv_',nn,3,sep='')]<-sum(tt3)
    cc
  }
}