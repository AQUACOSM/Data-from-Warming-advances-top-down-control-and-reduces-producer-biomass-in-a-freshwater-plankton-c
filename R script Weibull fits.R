library("cardidates", lib.loc="~/R/win-library/3.0")
library("boot", lib.loc="C:/Program Files/R/R-3.0.0/library")

#note: area under the curve analysis are done in the excelfiles with the modelfits from these respective weibull analysis

## Rotifers ----
data=read.csv('weibull_data_rotifers.csv')

#create dataframe for cardinal date output of weibull analysis
lLIMNO=c('X1', 'X2', 'X3', 'X5', 'X6', 'X7', 'X8', 'X9') 
lFITS=c('treatment','peakid', 'end1', 't1Mid','t1Begin','t1End')
cardates = data.frame(matrix(NA,length(lLIMNO),length(lFITS)))
colnames(cardates) = lFITS
rownames(cardates) = lLIMNO
cardates$treatment=c('warm','control', 'control', 'warm', 'warm', 'control', 'warm', 'control')

#create dataframe for the modelfits of the weibull analysis (for area under the curve calculations in excel)
modelfits=data.frame(matrix(NA, nrow=nrow(data), ncol=length(lLIMNO)))
xwaarden=data.frame(matrix(NA, nrow=nrow(data), ncol=length(lLIMNO)))
colnames(modelfits)=lLIMNO
colnames(xwaarden)=c('x-L1', 'x-L2', 'x-L3', 'x-L5', 'x-L6', 'x-L7', 'x-L8', 'x-L9')

#run the weibull analysis
nRUN=0
for(sLIMNO in lLIMNO){
  nRUN=nRUN+1
  peaks=peakwindow(data$date, data[,which(colnames(data)==sLIMNO,)])
  plot(peaks, main=sLIMNO)
  id=max(unique(peaks$peakid))
  cardates[nRUN,c(2)] = id
  enddate=data[which(rownames(data)==max(peaks$smd.indices)),1]
  cardates[nRUN,c(3)] = enddate
  set=data
  res=fitweibull6(set$date, set[,which(colnames(set)==sLIMNO)],linint = 0)
  card=CDW(res)$x
  cardates[nRUN,c(4:6)] = card 
  plot(res)
  modelfits[,which(colnames(modelfits)==sLIMNO)]=(res$fit$f*res$ymax)
  xwaarden[,which(colnames(modelfits)==sLIMNO)]=res$fit$x
}

#calculate differences in cardinal dates with t.tests
warm<-subset(cardates, treatment=="warm") 
control<-subset(cardates, treatment=='control')
t.test(control$t1Begin, warm$t1Begin)
t.test(control$t1Mid, warm$t1Mid)
t.test(control$t1End, warm$t1End)

#calculate the shift in days for the cardinal date
mean(control$t1End)-mean(warm$t1End)

## Cladocera ------------------------------------------------------------------------------
data=read.csv('weibull_data_cladocera.csv')

#create dataframe for cardinal date output of weibull analysis
lLIMNO=c('X1', 'X2', 'X3', 'X5', 'X6', 'X7', 'X8', 'X9')
lFITS=c('treatment','peakid', 'end1', 't1Mid','t1Begin','t1End')
cardates = data.frame(matrix(NA,length(lLIMNO),length(lFITS)))
colnames(cardates) = lFITS
rownames(cardates) = lLIMNO
cardates$treatment=c('warm','control', 'control', 'warm', 'warm', 'control', 'warm', 'control')

#create dataframe for the modelfits of the weibull analysis (for area under the curve calculations in excel)
modelfits=data.frame(matrix(NA, nrow=nrow(data), ncol=length(lLIMNO)))
xwaarden=data.frame(matrix(NA, nrow=nrow(data), ncol=length(lLIMNO)))
colnames(modelfits)=lLIMNO
colnames(xwaarden)=c('x-L1', 'x-L2', 'x-L3', 'x-L5', 'x-L6', 'x-L7', 'x-L8', 'x-L9')

#run the weibull analysis
nRUN=0
for(sLIMNO in lLIMNO){
  nRUN=nRUN+1
  peaks=peakwindow(data$date, data[,which(colnames(data)==sLIMNO,)])
  plot(peaks, main=sLIMNO)
  id=max(unique(peaks$peakid))
  cardates[nRUN,c(2)] = id
  enddate=data[which(rownames(data)==max(peaks$smd.indices)),1]
  cardates[nRUN,c(3)] = enddate
  set=data
  res=fitweibull6(set$date, set[,which(colnames(set)==sLIMNO)],linint = 0)
  card=CDW(res)$x
  cardates[nRUN,c(4:6)] = card 
  plot(res)
  modelfits[,which(colnames(modelfits)==sLIMNO)]=(res$fit$f*res$ymax)
  xwaarden[,which(colnames(modelfits)==sLIMNO)]=res$fit$x
}

#determine whether cardinal dates differ between treatments
warm<-subset(cardates, treatment=="warm") 
control<-subset(cardates, treatment=='control')
t.test(control$t1Begin, warm$t1Begin)
t.test(control$t1Mid, warm$t1Mid)
t.test(control$t1End, warm$t1End)

#calculate difference in cardinal dates
mean(control$t1Begin)-mean(warm$t1Begin)
mean(control$t1Mid)-mean(warm$t1Mid)

## Copepods ------------------------------------------------------------------------------
data=read.csv('weibull_data_copepods.csv')

#create dataframe for cardinal date output of weibull analysis
lLIMNO=c('X1', 'X2','X3', 'X5', 'X6', 'X7', 'X8', 'X9') 
lFITS=c('treatment','peakid', 'end1', 't1Mid','t1Begin','t1End')
cardates = data.frame(matrix(NA,length(lLIMNO),length(lFITS)))
colnames(cardates) = lFITS
rownames(cardates) = lLIMNO
cardates$treatment=c('warm','control', 'control', 'warm', 'warm', 'control', 'warm', 'control')

#create dataframe for the modelfits of the weibull analysis (for area under the curve calculations in excel)
modelfits=data.frame(matrix(NA, nrow=nrow(data), ncol=length(lLIMNO)))
xwaarden=data.frame(matrix(NA, nrow=nrow(data), ncol=length(lLIMNO)))
colnames(modelfits)=lLIMNO
colnames(xwaarden)=c('x-L1', 'x-L2', 'x-L3', 'x-L5', 'x-L6', 'x-L7', 'x-L8', 'x-L9')

#run the weibull analysis
nRUN=0
for(sLIMNO in lLIMNO){
  nRUN=nRUN+1
  peaks=peakwindow(data$date, data[,which(colnames(data)==sLIMNO,)])
  plot(peaks, main=sLIMNO)
  id=max(unique(peaks$peakid))
  cardates[nRUN,c(2)] = id
  enddate=data[which(rownames(data)==max(peaks$smd.indices)),1]
  cardates[nRUN,c(3)] = enddate
  set=data
  res=fitweibull6(data$date, data[,which(colnames(data)==sLIMNO)])
  card=CDW(res)$x
  cardates[nRUN,c(4:6)] = card 
  plot(res)
  modelfits[,which(colnames(modelfits)==sLIMNO)]=(res$fit$f*res$ymax)
  xwaarden[,which(colnames(modelfits)==sLIMNO)]=res$fit$x
}

#determine whether cardinal dates differ between treatments
warm<-subset(cardates, treatment=="warm") 
control<-subset(cardates, treatment=='control')
t.test(control$t1Begin, warm$t1Begin)
t.test(control$t1Mid, warm$t1Mid)
t.test(control$t1End, warm$t1End)

#calculate difference in cardinal dates
mean(control$t1Mid)-mean(warm$t1Mid)

## Chlorophyll-a (spring bloom) ------------------------------------------------------------------------------
data=read.csv('weibull_data_chla.csv')

#create dataframe for cardinal date output of weibull analysis
lLIMNO=c('X1', 'X2','X3', 'X5', 'X6', 'X7', 'X8', 'X9') 
lFITS=c('treatment','peakid', 'end1', 't1Mid','t1Begin','t1End')
cardates = data.frame(matrix(NA,length(lLIMNO),length(lFITS)))
colnames(cardates) = lFITS
rownames(cardates) = lLIMNO
cardates$treatment=c('warm','control', 'control', 'warm', 'warm', 'control', 'warm', 'control')

#create dataframe for the modelfits of the weibull analysis (for area under the curve calculations in excel)
modelfits=data.frame(matrix(NA, nrow=nrow(data), ncol=length(lLIMNO)))
xwaarden=data.frame(matrix(NA, nrow=nrow(data), ncol=length(lLIMNO)))
colnames(modelfits)=lLIMNO
colnames(xwaarden)=c('x-L1', 'x-L2', 'x-L3', 'x-L5', 'x-L6', 'x-L7', 'x-L8', 'x-L9')

#run the weibull analysis
nRUN=0
for(sLIMNO in lLIMNO){
  nRUN=nRUN+1
  peaks=peakwindow(data$date, data[,which(colnames(data)==sLIMNO,)])
  plot(peaks, main=sLIMNO)
  id=max(unique(peaks$peakid))
  cardates[nRUN,c(2)] = id
  enddate=data[which(rownames(data)==max(peaks$smd.indices)),1]
  cardates[nRUN,c(3)] = enddate
  set=data
  res=fitweibull6(set$date, set[,which(colnames(set)==sLIMNO)],linint = 0)
  card=CDW(res)$x
  cardates[nRUN,c(4:6)] = card 
  plot(res)
  modelfits[,which(colnames(modelfits)==sLIMNO)]=(res$fit$f*res$ymax)
  xwaarden[,which(colnames(modelfits)==sLIMNO)]=res$fit$x
}

#determine whether cardinal dates differ between treatments
warm<-subset(cardates, treatment=="warm") 
control<-subset(cardates, treatment=='control')
t.test(control$t1Begin, warm$t1Begin)
t.test(control$t1Mid, warm$t1Mid)
t.test(control$t1End, warm$t1End)

## Chytrid prevalence ------------------------------------------------------------------------------
data=read.csv('weibull_data_chytrids.csv')

#create dataframe for cardinal date output of weibull analysis
lLIMNO=c('X1', 'X2','X3', 'X5', 'X6', 'X7', 'X8', 'X9')
lFITS=c('treatment','peakid', 'end1', 't1Mid','t1Begin','t1End')
cardates = data.frame(matrix(NA,length(lLIMNO),length(lFITS)))
colnames(cardates) = lFITS
rownames(cardates) = lLIMNO
cardates$treatment=c('warm','control', 'control', 'warm', 'warm', 'control', 'warm', 'control')

#create dataframe for the modelfits of the weibull analysis (for area under the curve calculations in excel)
modelfits=data.frame(matrix(NA, nrow=nrow(data), ncol=length(lLIMNO)))
xwaarden=data.frame(matrix(NA, nrow=nrow(data), ncol=length(lLIMNO)))
colnames(modelfits)=lLIMNO
colnames(xwaarden)=c('x-L1', 'x-L2', 'x-L3', 'x-L5', 'x-L6', 'x-L7', 'x-L8', 'x-L9')

#run the weibull analysis
nRUN=0
for(sLIMNO in lLIMNO){
  nRUN=nRUN+1
  peaks=peakwindow(data$date, data[,which(colnames(data)==sLIMNO,)])
  plot(peaks, main=sLIMNO)
  id=max(unique(peaks$peakid))
  cardates[nRUN,c(2)] = id
  enddate=data[which(rownames(data)==max(peaks$smd.indices)),1]
  cardates[nRUN,c(3)] = enddate
  set=data
  res=fitweibull6(set$date, set[,which(colnames(set)==sLIMNO)],linint = 0)
  card=CDW(res)$x
  cardates[nRUN,c(4:6)] = card 
  plot(res)
  modelfits[,which(colnames(modelfits)==sLIMNO)]=(res$fit$f*res$ymax)
  xwaarden[,which(colnames(modelfits)==sLIMNO)]=res$fit$x
}

#determine whether cardinal dates differ between treatments
warm<-subset(cardates, treatment=="warm") 
control<-subset(cardates, treatment=='control')
t.test(control$t1Begin, warm$t1Begin)
t.test(control$t1Mid, warm$t1Mid)
t.test(control$t1End, warm$t1End)

#calculate difference in cardinal dates
mean(control$t1End)-mean(warm$t1End)
