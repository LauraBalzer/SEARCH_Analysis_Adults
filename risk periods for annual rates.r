# By Yea-Hung Chen

# risk period 1
ii$y1start<-ii$start
ii$y1end<-ii$end-365.25*2

# risk period 2
ii$y2start<-ii$end-365.25*2
ii$y2end<-ii$end-365.25

# risk period 3
ii$y3start<-ii$end-365.25
ii$y3end<-ii$end

# start date for adult
ii$adult_start<-as.Date('1900-01-01')
is.na(ii$adult_start)<-TRUE
ii$adult_start<-if_else(ii$age_0<15,
                        ii$y1start+(15-ii$age_0)*365.25,
                        ii$adult_start)
is.na(ii$adult_start)<-(ii$adult_start<ii$y1start)
is.na(ii$adult_start)<-(ii$y3end<ii$adult_start)

# start date for entering risk
ii$enter_start<-pmax(ii$inmigrant_date_3,ii$adult_start,na.rm=TRUE)
ii$enter_start<-with(ii,if_else(is.na(enter_start),y1start-365.25*100,
                                enter_start))

# modify risk periods for period-1 entries
ii$en.iny1<-with(ii,(y1start<enter_start)&(enter_start<y1end))
ii$y1start[ii$en.iny1]<-ii$enter_start[ii$en.iny1]

# modify risk periods for period-2 entries
ii$en.iny2<-with(ii,(y2start<enter_start)&(enter_start<y2end))
ii$y2start[ii$en.iny2]<-ii$enter_start[ii$en.iny2]
is.na(ii$y1start)<-ii$en.iny2
is.na(ii$y1end)<-ii$en.iny2

# modify risk periods for period-3 entries
ii$en.iny3<-with(ii,(y3start<enter_start)&(enter_start<y3end))
ii$y3start[ii$en.iny3]<-ii$enter_start[ii$en.iny3]
is.na(ii$y1start)<-ii$en.iny3
is.na(ii$y2start)<-ii$en.iny3
is.na(ii$y1end)<-ii$en.iny3
is.na(ii$y2end)<-ii$en.iny3