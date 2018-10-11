# By Yea-Hung Chen

# function for defining survival variables
define_survival_variables<-function(dd=ii){
  # convert all end dates to days
  end<-grep('\\.end$|^end$',names(dd),value=TRUE)
  dd[,end]<-lapply(end,function(ee){
    if(class(dd[,ee])=='Date'){
      ee<-as.numeric(dd[,ee]-dd$start)
      is.na(ee)<-(ee<0) 
      ee
    }
  })
  # define end variables
  oend<-paste(ovar,'.end',sep='')
  cend<-c(paste(cvar,'.end',sep=''),'end')
  # define custom function for finding minimum
  cmin<-function(vv){
    if(sum(is.na(vv))==length(vv)){
      mm<-9999
      is.na(mm)<-TRUE
      mm
    } else { min(vv,na.rm=TRUE) }
  }
  # define survival variables
  dd$oend<-apply(dd[,oend,drop=FALSE],1,cmin)/365.25
  dd$cend<-apply(dd[,cend,drop=FALSE],1,cmin)/365.25
  dd$rend<-pmin(dd$oend,dd$cend,na.rm=TRUE)
  dd$o<-!is.na(dd$oend)&(dd$oend==dd$rend)
  # output data
  dd
}