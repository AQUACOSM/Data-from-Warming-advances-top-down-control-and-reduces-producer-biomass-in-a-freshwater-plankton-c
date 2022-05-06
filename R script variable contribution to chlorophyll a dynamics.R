library(car)
library(zoo)
library(MASS)
library(forecast)
library(lmerTest)
library(lme4)
library(LMERConvenienceFunctions)
library(arm)
library(stats)
library(FactoMineR)
library(relaimpo)
library(perm)
library(reshape2)

data.frame <- read.csv("data_variable contribution to chla.csv")

## produce ts series --------------

# get names of columns 
all.col.names <- colnames(data.frame)
#drop the first five 
all.col.names <- all.col.names[6:length(all.col.names)] 

#start and end might vary with other timeseries, check with funcion ts
ts.start <- c(2014,10)
ts.end <- c(2014,33)
limnotrons=c(1,2,3,5,6,7,8,9)

output=data.frame(matrix(NA, nrow = 8, ncol = length(colnames(data.frame))-2))
colnames(output)[length(colnames(data.frame))-4]='R2 full'
colnames(output)[length(colnames(data.frame))-3]='R2 final'
colnames(output)[length(colnames(data.frame))-2]='limnotron'
colnames(output)[length(colnames(data.frame))-5]='#variables'
treatment=c('control', 'warm')

# AIC selection for which variables to include in the analysis
for(n in limnotrons){
  df=subset(data.frame, data.frame$limnotron==n)
  
  # determine if splining is necessary and do it
  spline.cols=colnames(df)[colSums(is.na(df)) > 0]
  print(spline.cols)
  for(col.name in spline.cols) { 
    df[col.name]<-na.spline(zoo(df[col.name],1:24),na.rm=FALSE,maxgap=Inf)
  }
  
  #create timeseries of selected variable
  print(n)
  for(col.name in all.col.names) { 
    x <- ts(df[col.name], 
            start=ts.start, 
            end=ts.end, 
            frequency=52); 
    # create, or replace, a global variable using assign() 
    assign(paste0(col.name,".ts"),x, envir=.GlobalEnv); 
  } 
  
  #start with full model (self defined) and do a backward selection on AIC criteria
  full <- lm(chla.ts ~  DIN.ts + DIP.ts +  temperature.ts + Si.ts + prevalence.ts + zooplankton.ts, data=df)
  step <- stepAIC(full)

  #substract final variables from step AIC selection and write to output table
  y=rownames(summary(step)$coefficients)
  print(y[2:length(y)])
  var.sel=y[2:(length(colnames(output))+1)]
  nvar=which(limnotrons==n)
  output[nvar,]=var.sel
  output[nvar,length(colnames(data.frame))-5]=length(y)-1
  sFORM_EXPL=""
  for (b in 2:length(y)){
    sFORM_EXPL=paste(sFORM_EXPL,(y[b]), "+",sep=" ")
  }
  sFORM_EXPL=substr(sFORM_EXPL, 1, nchar(sFORM_EXPL)-1) 
  final=lm(as.formula(paste("chla~",sFORM_EXPL)), data=df )
  output[nvar,length(colnames(data.frame))-4]=summary(full)$adj.r.squared
  output[nvar,length(colnames(data.frame))-3]=summary(final)$adj.r.squared
  output[nvar, length(colnames(data.frame))-2]=n
}
rownames(output)=limnotrons

# testing of residuals (for autocorrelation) to see if ARIMA transformation is necessary
arima.parameters=data.frame(matrix(NA, nrow = 8, ncol=4))
colnames(arima.parameters)=c('limnotron', 'p', 'd', 'q')
arima.parameters[,1]= limnotrons

for(n in limnotrons){
  df=subset(data.frame, data.frame$limnotron==n)
  # determine if splining is necessary and do it
  spline.cols=colnames(df)[colSums(is.na(df)) > 0]
  print(spline.cols)
  for(col.name in spline.cols) { 
    df[col.name]<-na.spline(zoo(df[col.name],1:24),na.rm=FALSE,maxgap=Inf)
  }
  print(n)
  for(col.name in all.col.names) { 
    x <- ts(df[col.name], 
            start=ts.start, 
            end=ts.end, 
            frequency=52); 
    # create, or replace, a global variable using assign() 
    assign(paste0(col.name,".ts"),x, envir=.GlobalEnv); 
  } 
  selectedmodel <- lm(chla.ts ~  DIN.ts + DIP.ts +  temperature.ts + Si.ts + prevalence.ts + zooplankton.ts, data=df)
  #print(summary(full))
  res <- residuals.lm(selectedmodel)
  plot(res, main=n)
  coeff <- coefficients(selectedmodel)
  tsres <- ts((res), start = ts.start, frequency = 52) 
  plot(tsres, main=n)
  acf(tsres)
  z=auto.arima(res)
  print(n)
  print(z)
  nvar=which(limnotrons==n)
  arima.parameters[nvar,2]=z$arma[1]
  arima.parameters[nvar,3]=z$arma[6]
  arima.parameters[nvar,4]=z$arma[2]
}
# based on this output ARIMA transform limnotron 1 and 5.

data.arma=data.frame(matrix(NA, nrow= length(data.frame[,1]), ncol=length(data.frame)))
colnames(data.arma)=colnames(data.frame)
data.arma$limnotron=data.frame$limnotron
data.arma$date=data.frame$date
data.arma$week=data.frame$week
data.arma$day=data.frame$day
data.arma$treatment=data.frame$treatment

## ARIMA transformation ------------------------
#produces graphs of all time series (main is limnotron number). In red are ARIMA transformed timeseries (when necessary)
for (n in limnotrons){  
  df=subset(data.frame, limnotron==n)
  spline.cols=colnames(df)[colSums(is.na(df)) > 0]
  for(col.name in spline.cols) { 
    df[col.name]<-na.spline(zoo(df[col.name],1:24),na.rm=FALSE,maxgap=Inf)
  }
  if(sum(arima.parameters[which(arima.parameters[,1]==n),])-n>0){
    p=subset(arima.parameters, limnotron==n)[,2]
    d=subset(arima.parameters, limnotron==n)[,3]
    q=subset(arima.parameters, limnotron==n)[,4]
    for (col.name in all.col.names){
      x <- ts(df[col.name], 
              start=ts.start, 
              end=ts.end, 
              frequency=52); # measured on a weekly basis, so 52 
      z=arima(x, order = c (p,d,q)) 
      x_fit<-fitted.values(z)
      plot(x, main=n)
      lines(x_fit, col="red")
      assign(paste0(col.name, '.ts'),x, envir=.GlobalEnv)
      data.arma[which(data.arma$limnotron==n),which(colnames(data.arma)==col.name)]=x 
    } 
  }else{
    for (col.name in all.col.names){
      x <- ts(df[col.name], 
              start=ts.start, 
              end=ts.end, 
              frequency=52); 
      plot(x, main=n)
      assign(paste0(col.name, '.ts'),x, envir=.GlobalEnv)
      data.arma[which(data.arma$limnotron==n),which(colnames(data.arma)==col.name)]=x 
    }
  }
}

## run the linear model with the ARIMA transformed timeserie data ---------------------------

#create dataframe for final r2
r2.selected=data.frame(matrix(NA, nrow=8, ncol=1))
colnames(r2.selected)='r2.selected'
rownames(r2.selected)=limnotrons

for (n in limnotrons){
  df=subset(data.arma, data.arma$limnotron==n)
  sel.model=lm(chla ~  DIN + DIP +  temperature + Si + prevalence + zooplankton, data=df)
  print(n)
  print(summary(sel.model))
  r2.selected[which(rownames(r2.selected)==n),1]=summary(sel.model)$adj.r.squared
}

#create data.frame for coefficients from model (per variable)
coef.mod.sel=data.frame(matrix(NA, nrow=8, ncol=9))
colnames(coef.mod.sel)=c('limnotron', 'treatment', '(Intercept)', 'DIN', 'DIP', 'Si', 'temperature', 'prevalence', 'zooplankton' )
coef.mod.sel[,1]=limnotrons
treatment=c('warm', 'control', 'control', 'warm', 'warm', 'control', 'warm', 'control')
coef.mod.sel[,2]=treatment

for(n in limnotrons){
  df=subset(data.arma, data.arma$limnotron==n)
  full <- lm(chla ~  DIN+ DIP + Si+  temperature  + prevalence + zooplankton, data=df)
  coef.data=matrix((coef(full)))
  rownames(coef.data)=names(coef(full))
  for (x in rownames(coef.data)){
    coef.mod.sel[which(coef.mod.sel[,1]==n),which(colnames(coef.mod.sel)==x)]=coef.data[which(rownames(coef.data)==x),]
  }
}

#test if coefficients differ between treatments
coef.mod.sel[9,1]='ttest'
coef.mod.sel[10,1]='average control'
coef.mod.sel[11,1]='stdev control'
coef.mod.sel[12,1]='average warm'
coef.mod.sel[13,1]='stdev warm'
warm=subset(coef.mod.sel, treatment=='warm')
control=subset(coef.mod.sel, treatment=='control')
for(col in colnames(coef.mod.sel[3:length(coef.mod.sel)])){
  ttest=t.test(warm[col], control[col])
  coef.mod.sel[9,col]=ttest$p.value
  coef.mod.sel[10,col]=mean(control[1:4, col])
  coef.mod.sel[11,col]=sd(control[1:4, col])
  coef.mod.sel[12,col]=mean(warm[1:4, col])
  coef.mod.sel[13,col]=sd(warm[1:4, col])
}

#determinde contribution variabeles to adj.r2

#again create a data.frame for output
lHEAD_OUT=colnames(coef.mod.sel)[4:length(colnames(coef.mod.sel))]
dfOUT=as.data.frame(matrix(NA,length(unique(data.arma$limnotron)),length(lHEAD_OUT)))
colnames(dfOUT) = lHEAD_OUT
rownames(dfOUT) = limnotrons

#control vars
for(nLIMNO in limnotrons){
  vMODEL_SEL1 =lHEAD_OUT
  if(sum(is.na(vMODEL_SEL1))==0){
    vMODEL_SEL = vMODEL_SEL1
  }else{
    vMODEL_SEL = vMODEL_SEL1[-which(is.na(vMODEL_SEL1))]
  }

  
  #build a formula for model
  #---extract model parameters
  sFORM_R2=""
  for (sVAR in vMODEL_SEL){
    if(is.na(sVAR)==FALSE){
      sFORM_R2=paste(sFORM_R2,sVAR, "+",sep=" ")
    }
  }
  sFORM_R2=substr(sFORM_R2, 1, nchar(sFORM_R2)-1) 
  formR2 = as.formula(paste('chla', ' ~ ', sFORM_R2, sep=""))
  
  #now run the model
  lmR2_LIMNO = lm(formR2, data=data.arma[which(data.arma$limnotron == nLIMNO),])
  #calculate relative importance
  if(length(vMODEL_SEL) >1){
    relimpLIMNO = calc.relimp.lm(lmR2_LIMNO, rela=FALSE)
    dfRELIMP = as.data.frame(relimpLIMNO$lmg)
  }else{
    dfRELIMP = as.data.frame(summary(lmR2_LIMNO)$adj.r.squared) 
    rownames(dfRELIMP)=vMODEL_SEL
  }
  for(nROW in 1:nrow(dfRELIMP)){
    dfOUT[which(rownames(dfOUT) %in% as.character(nLIMNO)),which(colnames(dfOUT) %in% rownames(dfRELIMP)[nROW])]=dfRELIMP[nROW,]
  }

}
dfOUT_ALL = cbind.data.frame(limnotrons,treatment, output[,c((length(colnames(data.arma))-5),(length(colnames(data.arma))-4),(length(colnames(data.arma))-3))], r2.selected,dfOUT)

#running the anova and create boxplot graphs for the results
anova_data=melt (dfOUT_ALL, id=c('limnotrons', 'treatment', '#variables', 'R2 full', 'R2 final', 'r2.selected'))
par(mar=c(12,5,1,1))
boxplot(value~variable*treatment, data=anova_data, las = 2)
an_dfOUT_ALL=aov(anova_data$value~anova_data$variable*anova_data$treatment)
print(summary(an_dfOUT_ALL))
posthoc=TukeyHSD(x=an_dfOUT_ALL, 'anova_data$variable')
posthoc_data=as.data.frame(posthoc$'anova_data$variable')
posthocselect= posthoc_data[which(posthoc_data$'p adj'<0.05),]
boxplot(value~variable, data=anova_data, las = 2)
control=subset(dfOUT_ALL, treatment=='control')
warm=subset(dfOUT_ALL, treatment=='warm')

#create dataframe which summarizes the output (=Table 3 in manuscript)
treatmentcontribution=as.data.frame(matrix(NA,nrow=7, ncol=4))
colnames(treatmentcontribution)=c('mean control', 'sd control', 'mean warm', 'sd warm')
rownames(treatmentcontribution)=colnames(dfOUT_ALL)[6:12]
for (variable in rownames(treatmentcontribution)){
  treatmentcontribution[variable, 1]=mean(control[1:4,variable])
  treatmentcontribution[variable, 2]=sd(control[1:4,variable])
  treatmentcontribution[variable, 3]=mean(warm[1:4,variable])
  treatmentcontribution[variable, 4]=sd(warm[1:4,variable])
}
rownames(treatmentcontribution)[1]='Total model'
print(treatmentcontribution)

