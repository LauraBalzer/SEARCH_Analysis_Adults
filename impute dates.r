# By Yea-Hung Chen

# function for imputing dates
impute_dates<-function(end){
  # mark missing dates
  ii[,paste(end,'.na',sep='')]<-lapply(end,function(ee){
    vv<-gsub('\\.end','',ee)
    is.na(ii[,ee])&ii[,vv]&!is.na(ii[,vv])})
  # impute dates
  ii[,end]<-lapply(end,function(ee){
    if_else(ii[,paste(ee,'.na',sep='')],
            ii$start+(ii$end-ii$start)/2,
            ii[,ee])
  })
  # output data
  ii
}