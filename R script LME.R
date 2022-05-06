library("stats", lib.loc="C:/Program Files/R/R-3.0.0/library")
library(afex)
library(nlme)
options(contrasts=c("contr.sum","contr.poly"))

### Chla (phytoplankton abundance) ---------
data=read.csv('biweekly_data.csv')
data[16:26] <- list(NULL)

#create log and sqrt transformations to be tested
logchla=log(data$chla+.0001)
sqrtchla=sqrt(data$chla)
data=cbind(data, logchla, sqrtchla)

#create dataframe for testing for normality
lDATES =  c(as.character(unique(data$date)))
lHEADER = c('chla', 'logchla', 'sqrtchla')
NORM.abundance = data.frame(matrix(0,length(lDATES),length(lHEADER)))
colnames(NORM.abundance) = lHEADER
rownames(NORM.abundance) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, na.rm=TRUE, date==sDATE)
  if(var(x$chla, na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$chla)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,1] = nSW_P
  }
  if(var(x$logchla, na.rm=TRUE) !=0){  
    nSW_P = shapiro.test(x$logchla)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrtchla, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrtchla)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,3] = nSW_P
  }
}
colSums(NORM.abundance<0.05)

#create dataframe for testing homegeneity of variance
abundance.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(abundance.VARTEST) = lHEADER
rownames(abundance.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='cold')

#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  
  nSW_P = var.test(x$chla, y$chla, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,1] = nSW_P  
  
  nSW_P = var.test(x$logchla, y$logchla, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$sqrtchla, y$sqrtchla, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,3] = nSW_P 
}
colSums(abundance.VARTEST<0.05) 

###based on these tests use logtranformed data 

#create dataframe
VAR <- data$logchla
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)

#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("chlorophyll")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

### Phytoplankton size classes ----------
data=read.csv('biweekly_data.csv')
data[18:26] <- list(NULL)

#sum prokaryote and eukaryotic phytoplankton groups for each sizeclass to get overall abundance per sizeclass
data[,18]=data$eu02+data$pro02
colnames(data)[18]='02'
data[,19]=data$eu230+data$pro230
colnames(data)[19]='230'
data[,20]=data$eu3085+data$pro3085
colnames(data)[20]='3085'

#create log and sqrt transformations to be tested
log02=log(data$'02'+.0001)
log230=log(data$'230'+.0001)
log3085=log(data$'3085'+.0001)
sqrt02=sqrt(data$'02')
sqrt230=sqrt(data$'230')
sqrt3085=sqrt(data$'3085')
data=cbind(data, sqrt02, sqrt230, sqrt3085, log02, log230, log3085)

#create dataframe for testing for normality
lDATES =  c(as.character(unique(data$date)))
lHEADER = c('02', 'log02', 'sqrt02', '230', 'log230', 'sqrt230', '3085', 'log3085','sqrt3085')
NORM.abundance = data.frame(matrix(0,length(lDATES),length(lHEADER)))
colnames(NORM.abundance) = lHEADER
rownames(NORM.abundance) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, na.rm=TRUE, date==sDATE)
  if(var(x$'02', na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$'02')$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,1] = nSW_P
  }
  if(var(x$log02, na.rm=TRUE) !=0){  
    nSW_P = shapiro.test(x$log02)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrt02, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrt02)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,3] = nSW_P
  }
  if(var(x$'230', na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$'230')$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,4] = nSW_P
  }
  if(var(x$log230, na.rm=TRUE) !=0){  
    nSW_P = shapiro.test(x$log230)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,5] = nSW_P
  }
  if(var(x$sqrt230, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrt230)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,6] = nSW_P
  }
  if(var(x$'3085', na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$'3085')$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,7] = nSW_P
  }
  if(var(x$log3085, na.rm=TRUE) !=0){  
    nSW_P = shapiro.test(x$log3085)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,8] = nSW_P
  }
  if(var(x$sqrt3085, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrt3085)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,9] = nSW_P
  }
}
colSums(NORM.abundance<0.05)

#create dataframe for testing homegeneity of variance
abundance.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(abundance.VARTEST) = lHEADER
rownames(abundance.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='cold')

#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  
  nSW_P = var.test(x$'02', y$'02', na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,1] = nSW_P  
  
  nSW_P = var.test(x$log02, y$log02, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$sqrt02, y$sqrt02, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,3] = nSW_P 
  
  nSW_P = var.test(x$'230', y$'230', na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,4] = nSW_P  
  
  nSW_P = var.test(x$log230, y$log230, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,5] = nSW_P 
  
  nSW_P = var.test(x$sqrt230, y$sqrt230, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,6] = nSW_P 
  
  nSW_P = var.test(x$'3085', y$'3085', na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,7] = nSW_P  
  
  nSW_P = var.test(x$log3085, y$log3085, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,8] = nSW_P 
  
  nSW_P = var.test(x$sqrt3085, y$sqrt3085, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,9] = nSW_P 
}
colSums(abundance.VARTEST<0.05) 

###based on these tests use logtranformed data for all size classes

VAR <- data$log02
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("log02")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

VAR <- data$log230
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("log230")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

VAR <- data$log3085
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("log3085")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

### seston C:N:P Stoichiometry ---------
data=read.csv('weekly_data.csv')

#create log and sqrt transformations to be tested
logCN220=log(data$CN220+.0001)
sqrtCN220=sqrt(data$CN220)
logCP220=log(data$CP220+.0001)
sqrtCP220=sqrt(data$CP220)
logNP220=log(data$NP220+.0001)
sqrtNP220=sqrt(data$NP220)
data=cbind(data, logCN220, sqrtCN220, logCP220, logNP220, sqrtCP220, sqrtNP220)

#create dataframe for testing for normality
lDATES = c(unique(as.character(data$date)))
lHEADER = c('CN220', 'CP220', 'NP220', 'logCN220', 'logCP220', 'logNP220', 'sqrtCN220', 'sqrtCP220', 'sqrtNP220')
NORM.stoich = data.frame(matrix(0,length(lDATES),length(lHEADER)))
colnames(NORM.stoich) = lHEADER
rownames(NORM.stoich) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, date==sDATE)
  if(var(x$CN220)!=0){
    nSW_P = shapiro.test(x$CN220)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,1] = nSW_P
  }
  if(var(x$CP220) !=0){  
    nSW_P = shapiro.test(x$CP220)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,2] = nSW_P
  }
  if(var(x$NP220) !=0){
    nSW_P = shapiro.test(x$NP220)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,3] = nSW_P
  }
  if(var(x$logCN220) !=0){
    nSW_P = shapiro.test(x$logCN220)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,4] = nSW_P
  }
  if(var(x$logCP220) !=0){
    nSW_P = shapiro.test(x$logCP220)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,5] = nSW_P
  }
  if(var(x$logNP220) !=0){
    nSW_P = shapiro.test(x$logNP220)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,6] = nSW_P
  }
  if(var(x$sqrtCN220) !=0){
    nSW_P = shapiro.test(x$sqrtCN220)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,7] = nSW_P
  }
  if(var(x$sqrtCP220) !=0){
    nSW_P = shapiro.test(x$sqrtCP220)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,8] = nSW_P
  }  
  if(var(x$sqrtNP220) !=0){
    nSW_P = shapiro.test(x$sqrtNP220)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,9] = nSW_P
  }
}
colSums(NORM.stoich<0.05)

#create dataframe for testing homegeneity of variance
stoich.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(stoich.VARTEST) = lHEADER
rownames(stoich.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='control')

#variance test per timepoint for all transfomrations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  
  nSW_P = var.test(x$CN220, y$CN220)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  stoich.VARTEST[nDATE_LOC,1] = nSW_P  

  nSW_P = var.test(x$CP220, y$CP220)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  stoich.VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$NP220, y$NP220)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  stoich.VARTEST[nDATE_LOC,3] = nSW_P 
 
  nSW_P = var.test(x$logCN220, y$logCN220)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  stoich.VARTEST[nDATE_LOC,4] = nSW_P 
  
  nSW_P = var.test(x$logCP220, y$logCP220)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  stoich.VARTEST[nDATE_LOC,5] = nSW_P 
  
  nSW_P = var.test(x$logNP220, y$logNP220)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  stoich.VARTEST[nDATE_LOC,6] = nSW_P 
  
  nSW_P = var.test(x$sqrtCN220, y$sqrtCN220)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  stoich.VARTEST[nDATE_LOC,7] = nSW_P 
  
  nSW_P = var.test(x$sqrtCP220, y$sqrtCP220)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  stoich.VARTEST[nDATE_LOC,8] = nSW_P 
  
  nSW_P = var.test(x$sqrtNP220, y$sqrtNP220)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  stoich.VARTEST[nDATE_LOC,9] = nSW_P 
}
colSums(stoich.VARTEST<0.05)

### based on these tests, continue with raw data CN, log transformed CP and sqrt transformed NP

VAR <- data$CN220
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("CN")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

VAR <- data$logCP220
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("logCP")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

VAR <- data$sqrtNP220
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("sqrtNP")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

### DIN and DIP -------
data=read.csv('biweekly_data.csv')
data[16:25] <- list(NULL)


#create data transformations to be tested
logN=log(data$DIN+.0001)
logP=log(data$DIP+.0001)
sqrtN=sqrt(data$DIN)
sqrtP=sqrt(data$DIP)
frthrtN=data$DIN^0.25
frthrtP=data$DIP^0.25
data=cbind(data, logN, logP, sqrtN, sqrtP, frthrtN, frthrtP)

#create dataframe for testing for normality
lDATES = c(unique(as.character(data$date)))
lHEADER = c('DIN', 'DIP', 'logN', 'logP', 'sqrtN', 'sqrtP', 'frthrtN', 'frthrtP')
NORM.nuts = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(NORM.nuts) = lHEADER
rownames(NORM.nuts) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, date==sDATE)
  if(var(x$DIN)!=0){
    nSW_P = shapiro.test(x$DIN)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.nuts[nDATE_LOC,1] = nSW_P
  }
  if(var(x$DIP) !=0){  
    nSW_P = shapiro.test(x$DIP)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.nuts[nDATE_LOC,2] = nSW_P
  }
  if(var(x$logN) !=0){
    nSW_P = shapiro.test(x$logN)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.nuts[nDATE_LOC,3] = nSW_P
  }
  if(var(x$logP) !=0){
    nSW_P = shapiro.test(x$logP)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.nuts[nDATE_LOC,4] = nSW_P
  }
  if(var(x$sqrtN) !=0){
    nSW_P = shapiro.test(x$sqrtN)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.nuts[nDATE_LOC,5] = nSW_P
  }
  if(var(x$sqrtP) !=0){
    nSW_P = shapiro.test(x$sqrtP)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.nuts[nDATE_LOC,6] = nSW_P
  }
  if(var(x$frthrtN) !=0){
    nSW_P = shapiro.test(x$frthrtN)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.nuts[nDATE_LOC,7] = nSW_P
  }
  if(var(x$frthrtP) !=0){
    nSW_P = shapiro.test(x$frthrtP)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.nuts[nDATE_LOC,8] = nSW_P
  }
}
colSums(NORM.nuts<0.05, na.rm=TRUE) 

#create dataframe for testing homegeneity of variance
VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(VARTEST) = lHEADER
rownames(VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='cold')

#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  
  nSW_P = var.test(x$DIN, y$DIN, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  VARTEST[nDATE_LOC,1] = nSW_P  
  
  nSW_P = var.test(x$DIP, y$DIP, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$logN, y$logN, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  VARTEST[nDATE_LOC,3] = nSW_P 
  
  nSW_P = var.test(x$logP, y$logP, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  VARTEST[nDATE_LOC,4] = nSW_P   
  
  nSW_P = var.test(x$sqrtN, y$sqrtN, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  VARTEST[nDATE_LOC,5] = nSW_P   
  
  nSW_P = var.test(x$sqrtP, y$sqrtP, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  VARTEST[nDATE_LOC,6] = nSW_P   
  
  nSW_P = var.test(x$frthrtN, y$frthrtN, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  VARTEST[nDATE_LOC,7] = nSW_P 
  
  nSW_P = var.test(x$frthrtP, y$frthrtP, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  VARTEST[nDATE_LOC,8] = nSW_P
}
colSums(VARTEST<0.05, na.rm=TRUE)

## use raw data in analysis based on these tests

VAR <- data$DIN
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("DIN")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

VAR <- data$DIP
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("DIP")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

### Dissolved Si ----
data=read.csv('Si_data.csv')

#create log and sqrt transformations to be tested
logsi=log(data$Si+.0001)
sqrtsi=sqrt(data$Si)
data=cbind(data, logsi, sqrtsi)

#create dataframe for testing for normality
lDATES =  c(as.character(unique(data$date)))
lHEADER = c('Si', 'logSi', 'sqrtSi')
NORM.abundance = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(NORM.abundance) = lHEADER
rownames(NORM.abundance) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, na.rm=TRUE, date==sDATE)
  if(var(x$Si, na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$Si)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,1] = nSW_P
  }
  if(var(x$logsi, na.rm=TRUE) !=0){  
    nSW_P = shapiro.test(x$logsi)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrtsi, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrtsi)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,3] = nSW_P
  }
}
colSums(NORM.abundance<0.05) ## sqrt data

#create dataframe for testing homegeneity of variance
abundance.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(abundance.VARTEST) = lHEADER
rownames(abundance.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='control')

#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  
  nSW_P = var.test(x$Si, y$Si, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,1] = nSW_P  
  
  nSW_P = var.test(x$logsi, y$logsi, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$sqrtsi, y$sqrtsi, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,3] = nSW_P 
}
colSums(abundance.VARTEST<0.05)

## use sqrt data in analysis

VAR <- data$sqrtsi
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("Si")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

### Dissolved CO2 ---------
data=read.csv('weekly_data.csv')

#create data transformations to be tested
logc=log(data$pCO2+.0001) 
sqrtc=sqrt(data$pCO2)
frthrtc=(data$pCO2)^0.25
loglogc=log(data$pCO2+.0001)
data=cbind(data, logc, sqrtc, frthrtc, loglogc)

#create dataframe for testing for normality
lDATES =  c(as.character(unique(data$date)))
lHEADER = c('pCO2', 'logc', 'sqrtc', 'frthrtc', 'loglogc')
PVAL.tic = data.frame(matrix(0,length(lDATES),length(lHEADER)))
colnames(PVAL.tic) = lHEADER
rownames(PVAL.tic) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, date==sDATE)
  if(var(x$pCO2,na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$pCO2)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    PVAL.tic[nDATE_LOC,1] = nSW_P
  }
  if(var(x$logc) !=0){  
    nSW_P = shapiro.test(x$logc)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    PVAL.tic[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrtc) !=0){
    nSW_P = shapiro.test(x$sqrtc)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    PVAL.tic[nDATE_LOC,3] = nSW_P
  }
  if(var(x$frthrtc) !=0){
    nSW_P = shapiro.test(x$frthrtc)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    PVAL.tic[nDATE_LOC,4] = nSW_P
  }
  if(var(x$loglogc) !=0){
    nSW_P = shapiro.test(x$loglogc)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    PVAL.tic[nDATE_LOC,5] = nSW_P
  }
}
colSums(PVAL.tic<0.05) ## sqrt

#create dataframe for testing homegeneity of variance
TIC.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(TIC.VARTEST) = lHEADER
rownames(TIC.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='control')

#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  nSW_P = var.test(x$pCO2, y$pCO2)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  TIC.VARTEST[nDATE_LOC,1] = nSW_P  
  
  nSW_P = var.test(x$logc, y$logc)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  TIC.VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$sqrtc, y$sqrtc)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  TIC.VARTEST[nDATE_LOC,3] = nSW_P 
  
  nSW_P = var.test(x$frthrtc, y$frthrtc)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  TIC.VARTEST[nDATE_LOC,4] = nSW_P 
  
  nSW_P = var.test(x$loglogc, y$loglogc)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  TIC.VARTEST[nDATE_LOC,5] = nSW_P 
}
colSums(TIC.VARTEST<0.05) 

## based on these tests, use sqrt transformed data 

VAR <- data$sqrtc
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("sqrtc")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

### Light availability -----------
data=read.csv('weekly_data.csv')

#create log and sqrt transformations to be tested
logI50=log(data$I50+.0001) 
sqrtI50=sqrt(data$I50)
frthrtI50=(data$I50)^0.25
loglogI50=log(logI50+.0001)
data=cbind(data, logI50, sqrtI50, frthrtI50, loglogI50)

#create dataframe for testing for normality
lDATES =  c(as.character(unique(data$date)))
lHEADER = c('I50', 'logI50', 'sqrtI50', 'frthrtI50', 'loglogI50')
PVAL.I = data.frame(matrix(0,length(lDATES),length(lHEADER)))
colnames(PVAL.I) = lHEADER
rownames(PVAL.I) = lDATES

#create dataframe for testing homegeneity of variance
I.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(I.VARTEST) = lHEADER
rownames(I.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='control')

#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  nSW_P = var.test(x$I50, y$I50)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  I.VARTEST[nDATE_LOC,1] = nSW_P  
  
  nSW_P = var.test(x$logI50, y$logI50)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  I.VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$sqrtI50, y$sqrtI50)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  I.VARTEST[nDATE_LOC,3] = nSW_P 
  
  nSW_P = var.test(x$frthrtI50, y$frthrtI50)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  I.VARTEST[nDATE_LOC,4] = nSW_P 
  
  nSW_P = var.test(x$loglogI50, y$loglogI50)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  I.VARTEST[nDATE_LOC,5] = nSW_P 
}
colSums(I.VARTEST<0.05)

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, date==sDATE)
  if(var(x$I50)!=0){
    nSW_P = shapiro.test(x$I50)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    PVAL.I[nDATE_LOC,1] = nSW_P
  }
  if(var(x$logI50) !=0){  
    nSW_P = shapiro.test(x$logI50)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    PVAL.I[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrtI50) !=0){
    nSW_P = shapiro.test(x$sqrtI50)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    PVAL.I[nDATE_LOC,3] = nSW_P
  }
  if(var(x$frthrtI50) !=0){
    nSW_P = shapiro.test(x$frthrtI50)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    PVAL.I[nDATE_LOC,4] = nSW_P
  }
  if(var(x$loglogI50) !=0){
    nSW_P = shapiro.test(x$loglogI50)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    PVAL.I[nDATE_LOC,5] = nSW_P
  }
}
colSums(PVAL.I<0.05) 

## use raw data for analysis based on these tests

VAR <- data$I50
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("I50")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

### Zooplankton ----------
data=read.csv('LME_zooplankton_data.csv')

#create log and sqrt transformations to be tested
logrot=log(data$rotifers+.0001)
sqrtrot=sqrt(data$rotifers)
logclad=log(data$cladocerans+.0001)
sqrtclad=sqrt(data$cladocerans)
logcop=log(data$copepods+.0001)
sqrtcop=sqrt(data$copepods)
data=cbind(data, logrot, sqrtrot, logclad, sqrtclad, logcop, sqrtcop)

#create dataframe for testing for normality
lDATES =  c(as.character(unique(data$date)))
lHEADER = c('rotifers','logrot', 'sqrtrot', 'cladocerans', 'logclad', 'sqrtclad', 'copepods', 'logcop', 'sqrtcop')
NORM.stoich = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(NORM.stoich) = lHEADER
rownames(NORM.stoich) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, date==sDATE)
  if(var(x$rotifers, na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$rotifers)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,1] = nSW_P
  }
  if(var(x$logrot, na.rm=TRUE) !=0){  
    nSW_P = shapiro.test(x$logrot)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrtrot, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrtrot)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,3] = nSW_P
  }
  if(var(x$cladocerans, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$cladocerans)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,4] = nSW_P
  }
  if(var(x$logclad, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$logclad)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,5] = nSW_P
  }
  if(var(x$sqrtclad, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrtclad)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,6] = nSW_P
  }
  if(var(x$copepods, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$copepods)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,7] = nSW_P
  }
  if(var(x$logcop, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$logcop)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,8] = nSW_P
  }  
  if(var(x$sqrtcop, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrtcop)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.stoich[nDATE_LOC,9] = nSW_P
  }
}
colSums(NORM.stoich<0.05) 

#create dataframe for testing homegeneity of variance
stoich.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(stoich.VARTEST) = lHEADER
rownames(stoich.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='control')

#variance test per timepoint for all transfomrations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
      nSW_P = var.test(x$rotifers, y$rotifers)$p.value
      nDATE_LOC = which(sDATE == lDATES)
      stoich.VARTEST[nDATE_LOC,1] = nSW_P  

      nSW_P = var.test(x$logrot, y$logrot)$p.value
      nDATE_LOC = which(sDATE == lDATES)
      stoich.VARTEST[nDATE_LOC,2] = nSW_P 

      nSW_P = var.test(x$sqrtrot, y$sqrtrot)$p.value
      nDATE_LOC = which(sDATE == lDATES)
      stoich.VARTEST[nDATE_LOC,3] = nSW_P 

      nSW_P = var.test(x$cladocerans, y$cladocerans)$p.value
      nDATE_LOC = which(sDATE == lDATES)
      stoich.VARTEST[nDATE_LOC,4] = nSW_P 

      nSW_P = var.test(x$logclad, y$logclad)$p.value
      nDATE_LOC = which(sDATE == lDATES)
      stoich.VARTEST[nDATE_LOC,5] = nSW_P 

      nSW_P = var.test(x$sqrtclad, y$sqrtclad)$p.value
      nDATE_LOC = which(sDATE == lDATES)
      stoich.VARTEST[nDATE_LOC,6] = nSW_P 

      nSW_P = var.test(x$copepods, y$copepods)$p.value
      nDATE_LOC = which(sDATE == lDATES)
      stoich.VARTEST[nDATE_LOC,7] = nSW_P 

      nSW_P = var.test(x$logcop, y$logcop)$p.value
      nDATE_LOC = which(sDATE == lDATES)
      stoich.VARTEST[nDATE_LOC,8] = nSW_P 

      nSW_P = var.test(x$sqrtcop, y$sqrtcop)$p.value
      nDATE_LOC = which(sDATE == lDATES)
      stoich.VARTEST[nDATE_LOC,9] = nSW_P 

}
colSums(stoich.VARTEST<0.05, na.rm=TRUE) 

## based on these tests use log transformed data for rotifers and sqrt transformation for cladocerans and copepods

VAR <- data$logrot
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("logrot")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

VAR <- data$sqrtclad
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("logclad")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

VAR <- data$sqrtcop
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("logcop")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

### Chytrid prevalence -------
data=read.csv('biweekly_data.csv')
data[16:25] <- list(NULL)

#create log and sqrt transformations to be tested
logchy=log(data$chytrids+.0001)
sqrtchy=asin(data$chytrids/100)
data=cbind(data, logchy, sqrtchy)

#create dataframe for testing for normality
lDATES =  c(as.character(unique(data$date)))
lHEADER = c('chytrids', 'logchy', 'sqrtchy')
NORM.abundance = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(NORM.abundance) = lHEADER
rownames(NORM.abundance) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, na.rm=TRUE, date==sDATE)
  if(var(x$chytrids, na.rm=TRUE)!=0 & is.na(var(x$chytrids, na.rm=TRUE))==F){
    nSW_P = shapiro.test(x$chytrids)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,1] = nSW_P
  }
  if(var(x$logchy, na.rm=TRUE) !=0 & is.na(var(x$logchy, na.rm=TRUE))==F){  
    nSW_P = shapiro.test(x$logchy)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrtchy, na.rm=TRUE) !=0 & is.na(var(x$sqrtchy, na.rm=TRUE))==F){
    nSW_P = shapiro.test(x$sqrtchy)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,3] = nSW_P
  }
}
colSums(NORM.abundance<0.05, na.rm=TRUE)

#create dataframe for testing homegeneity of variance
abundance.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(abundance.VARTEST) = lHEADER
rownames(abundance.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='cold')


#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  if(sum(x$chytrids, na.rm=T)!=0 & sum(y$chytrids, na.rm=T)!=0){
    nSW_P = var.test(x$chytrids, y$chytrids, na.rm=TRUE)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    abundance.VARTEST[nDATE_LOC,1] = nSW_P
  }
  if(sum(x$logchy, na.rm=T)!=0 & sum(y$logchy, na.rm=T)!=0){
    nSW_P = var.test(x$logchy, y$logchy, na.rm=TRUE)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    abundance.VARTEST[nDATE_LOC,2] = nSW_P
  }
  if(sum(x$sqrtchy, na.rm=T)!=0 & sum(y$sqrtchy, na.rm=T)!=0){
    nSW_P = var.test(x$sqrtchy, y$sqrtchy, na.rm=TRUE)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    abundance.VARTEST[nDATE_LOC,3] = nSW_P 
  }
}
colSums(abundance.VARTEST<0.05, na.rm=TRUE) 

### based on these tests, use raw data in analysis

VAR <- data$chytrids
TREATMENT <- factor(data$treatment)
TIME <- data$date
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("chytrids")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

### Periphyton ----
data=read.csv('periphyton_data.csv')

#create log and sqrt transformations to be tested
log=log(data$periphyton+.0001)
sqrt=sqrt(data$periphyton)
data=cbind(data, log, sqrt)

#create dataframe for testing for normality
lDATES = c(unique(as.character(data$date)))
lHEADER = c('periphyton', 'log', 'sqrt')
NORM.abundance = data.frame(matrix(0,length(lDATES),length(lHEADER)))
colnames(NORM.abundance) = lHEADER
rownames(NORM.abundance) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, na.rm=TRUE, date==sDATE)
  if(var(x$periphyton, na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$periphyton)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,1] = nSW_P
  }
  if(var(x$log, na.rm=TRUE) !=0){  
    nSW_P = shapiro.test(x$log)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrt, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrt)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,3] = nSW_P
  }
}
colSums(NORM.abundance<0.05)

#create dataframe for testing homegeneity of variance
abundance.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(abundance.VARTEST) = lHEADER
rownames(abundance.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='control')

#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  
  nSW_P = var.test(x$periphyton, y$periphyton, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,1] = nSW_P  
  
  nSW_P = var.test(x$log, y$log, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$sqrt, y$sqrt, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,3] = nSW_P 
}
colSums(abundance.VARTEST<0.05, na.rm=TRUE) 

# based on these tests use logtransformed data in the analysis

VAR <- data$log
TREATMENT <- factor(data$treatment)
SUBJECT <- data$limnotron
TIME <- factor(data$date)
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
baseline <- lme(VAR ~ 1, random = ~1 | SUBJECT/TIME, data = dfREPANOVA, method = "ML", na.action=na.omit)
treatment <- lme(VAR ~ TREATMENT, random = ~1 | SUBJECT/TIME, data = dfREPANOVA, method = "ML", na.action=na.omit)
time <- lme(VAR ~ TREATMENT+TIME, random = ~1 | SUBJECT/TIME, data = dfREPANOVA, method = "ML", na.action=na.omit)
treatmenttime <- lme(VAR ~ TREATMENT*TIME, random = ~1 | SUBJECT/TIME, data = dfREPANOVA, method = "ML", na.action=na.omit)
print("logperiphyton")
print(anova(baseline, treatment, time, treatmenttime))

### Benthic algae -----
data=read.csv('epipelon_data.csv')

#create log and sqrt transformations to be tested
log=log(data$epipelon+.0001)
sqrt=sqrt(data$epipelon)
data=cbind(data, log, sqrt)

#create dataframe for testing for normality
lDATES =c(unique(as.character(data$date)))
lHEADER = c('epipelon', 'log', 'sqrt')
NORM.abundance = data.frame(matrix(0,length(lDATES),length(lHEADER)))
colnames(NORM.abundance) = lHEADER
rownames(NORM.abundance) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, na.rm=TRUE, date==sDATE)
  if(var(x$epipelon, na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$epipelon)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,1] = nSW_P
  }
  if(var(x$log, na.rm=TRUE) !=0){  
    nSW_P = shapiro.test(x$log)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrt, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrt)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,3] = nSW_P
  }
}
colSums(NORM.abundance<0.05)

#create dataframe for testing homegeneity of variance
abundance.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(abundance.VARTEST) = lHEADER
rownames(abundance.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='control')

#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  
  nSW_P = var.test(x$epipelon, y$epipelon, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,1] = nSW_P  
  
  nSW_P = var.test(x$log, y$log, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$sqrt, y$sqrt, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,3] = nSW_P 
}
colSums(abundance.VARTEST<0.05, na.rm=TRUE)

# based on these tests use log transformed data in the analysis

VAR <- data$log
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#remove timepoint to work around LME error message
dfREPANOVA=dfREPANOVA[-c(which(dfREPANOVA$TIME=='17-3-2014')), ]

#run the actual analysis
baseline <- lme(VAR ~ 1, random = ~1 | SUBJECT/TIME, data = dfREPANOVA, method = "ML", na.action=na.omit)
treatment <- lme(VAR ~ TREATMENT, random = ~1 | SUBJECT/TIME, data = dfREPANOVA, method = "ML", na.action=na.omit)
time <- lme(VAR ~ TREATMENT+TIME, random = ~1 | SUBJECT/TIME, data = dfREPANOVA, method = "ML", na.action=na.omit)
treatmenttime <- lme(VAR ~ TREATMENT*TIME, random = ~1 | SUBJECT/TIME, data = dfREPANOVA, method = "ML", na.action=na.omit)
print("logepipelon")
print(anova(baseline, treatment, time, treatmenttime))

### Filamentous algae (PVI) ----
data=read.csv('Floating_and_filamentousalgae_data.csv')

#create log and sqrt transformations to be tested
logpvi=log(data$PVI+.0001)
sqrtpvi=sqrt(data$PVI)
data=cbind(data, logpvi, sqrtpvi)

#create dataframe for testing for normality
lDATES = c('27-jun',  '2-jul', '10-jul', '23-jul', '1-aug',  '8-aug', '19-sep')
lHEADER = c('PVI', 'logpvi', 'sqrtpvi')
NORM.abundance = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(NORM.abundance) = lHEADER
rownames(NORM.abundance) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, na.rm=TRUE, date==sDATE)
  if(var(x$PVI, na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$PVI)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,1] = nSW_P
  }
  if(var(x$logpvi, na.rm=TRUE) !=0){  
    nSW_P = shapiro.test(x$logpvi)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrtpvi, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrtpvi)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,3] = nSW_P
  }
}
colSums(NORM.abundance<0.05) ## sqrt data

#create dataframe for testing homegeneity of variance
abundance.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(abundance.VARTEST) = lHEADER
rownames(abundance.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='control')

#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  
  nSW_P = var.test(x$PVI, y$PVI, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,1] = nSW_P  
  
  nSW_P = var.test(x$logpvi, y$logpvi, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$sqrtpvi, y$sqrtpvi, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,3] = nSW_P 
}
colSums(abundance.VARTEST<0.05)

## based on these tests, use sqrt data in analysis

VAR <- data$sqrtpvi
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("sqrt filamentous algae")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))

### Floating algae (cover) ----
data=read.csv('Floating_and_filamentousalgae_data.csv')

#create log and sqrt transformations to be tested
logcover=log(data$cover+.0001)
sqrtcover=sqrt(data$cover)
data=cbind(data, logcover, sqrtcover)

#create dataframe for testing for normality
lDATES = c('27-jun',  '2-jul', '10-jul', '23-jul', '1-aug',  '8-aug', '19-sep')
lHEADER = c('cover', 'logcover', 'sqrtcover')
NORM.abundance = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(NORM.abundance) = lHEADER
rownames(NORM.abundance) = lDATES

#test of normality per timepoint 
for(sDATE in lDATES){
  x=subset(data, na.rm=TRUE, date==sDATE)
  if(var(x$cover, na.rm=TRUE)!=0){
    nSW_P = shapiro.test(x$cover)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,1] = nSW_P
  }
  if(var(x$logcover, na.rm=TRUE) !=0){  
    nSW_P = shapiro.test(x$logcover)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,2] = nSW_P
  }
  if(var(x$sqrtcover, na.rm=TRUE) !=0){
    nSW_P = shapiro.test(x$sqrtcover)$p.value
    nDATE_LOC = which(sDATE == lDATES)
    NORM.abundance[nDATE_LOC,3] = nSW_P
  }
}
colSums(NORM.abundance<0.05) ## sqrt data

#create dataframe for testing homegeneity of variance
abundance.VARTEST = data.frame(matrix(NA,length(lDATES),length(lHEADER)))
colnames(abundance.VARTEST) = lHEADER
rownames(abundance.VARTEST) = lDATES
warm<-subset(data, treatment=="warm") 
control<-subset(data, treatment=='control')

#variance test per timepoint for all transformations
for(sDATE in lDATES){
  x=subset(control, date==sDATE)
  y=subset(warm, date==sDATE)
  
  nSW_P = var.test(x$cover, y$cover, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,1] = nSW_P  
  
  nSW_P = var.test(x$logcover, y$logcover, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,2] = nSW_P 
  
  nSW_P = var.test(x$sqrtcover, y$sqrtcover, na.rm=TRUE)$p.value
  nDATE_LOC = which(sDATE == lDATES)
  abundance.VARTEST[nDATE_LOC,3] = nSW_P 
}
colSums(abundance.VARTEST<0.05)

## based on these results, use sqrt data in the analysis

VAR <- data$sqrtcover
TREATMENT <- factor(data$treatment)
TIME <- factor(data$date)
SUBJECT <- data$limnotron
dfREPANOVA <- data.frame(VAR, TREATMENT, TIME, SUBJECT)
#run the actual analysis
full=lme(VAR~TREATMENT*TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
interactions=lme(VAR~TREATMENT+TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
treatment=lme(VAR~TIME, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
time=lme(VAR~TREATMENT, random=~1|SUBJECT/TIME, data=dfREPANOVA, na.action=na.omit, method='ML')
print("sqrt floating algae")
print(anova(full, interactions, treatment))
print("for time effect")
print(anova(interactions, time))