multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

###------------ GUS ---==================
###--------------------================== Script to test whether foals from treated mares were affected =========--------
###----- 21.07.2020 ---==================
require(ggplot2)
require(Hmisc)
require(lme4)
require(lmerTest)
require(viridis)
require(fitdistrplus)
require(nlme)
require(geeM)

theme_set(theme_bw())
cols_tw = viridis_pal(option='D')(10)

setwd("~/Documents/INRA/HorseDrugResistance/IVMilk")

###----- Input files
# Faecal Egg Count
fec = read.csv(file="foal.csv",header=TRUE,sep=";")
colnames(fec) = c('Id','SamplingDate','AnalysisDate','eggsStr','FECs','eggsParascaris','FECp')
head(fec)
#     Id SamplingDate AnalysisDate eggsStr FECs eggsParascaris FECp
# 1 W650   28/05/2019   03/06/2019       0    0              0    0
# 2 W654   28/05/2019   03/06/2019       0    0              0    0
# 3 W655   28/05/2019   03/06/2019       0    0              0    0
# 4 W662   28/05/2019   03/06/2019       0    0              0    0
# 5 W675   28/05/2019   03/06/2019       8  400              0    0
# 6 W676   28/05/2019   03/06/2019       0    0              0    0
fec$Date= format(as.Date(fec$AnalysisDate, format ="%d/%m/%Y"),"%m/%d")

# Pasture
pasture = read.csv(file="group.csv",header=TRUE,sep=";")
colnames(pasture) = c('Id','Type','Pasture','Date','DateOfBirth')
head(pasture)
#      Id Type Pasture       Date DateOfBirth
# 1 MW585 FOAL      PA 25/05/2019  25/05/2019
# 2 MW586 FOAL      PA 19/06/2019  01/06/2019
# 3 MW587 FOAL      PA 04/06/2019  01/06/2019
# 4 MW588 FOAL      PA 08/06/2019  07/06/2019
# 5 MW589 FOAL      PA 24/06/2019  09/06/2019
# 6 MW590 FOAL      PN 19/06/2019  09/06/2019

fecpa = merge(fec,pasture[,c('Id', "Pasture", "Type", "DateOfBirth")],by='Id')

# Families
fam = read.csv(file = './families.csv', header=T, sep=';')
head(fam)
#   MARE STALLION DateOfBirth  FOAL SEX FoalWeight MareWeight
# 1 W728    MW511  09/05/2019  W756   F       31.2      349.2
# 2 W702    MW438  01/06/2019 MW587   M       21.8      263.0
# 3 W685    MW546  31/05/2019  W759   F       24.6      295.0
# 4 W676    MW546  25/05/2019 MW585   M       26.0      290.0
# 5 W707    MW438  01/06/2019 MW586   M       20.0      236.8
# 6 W678    MW438  28/05/2019  W757   F       26.0      271.0

#Discard W728 which was not part of the experiment
fam = fam[fam$MARE!='W728',]
fam$MARE = factor(fam$MARE)
fam$FOAL = factor(fam$FOAL)

### Convert variables
fecpa$AnalysisDate= format(as.Date(fecpa$AnalysisDate, format="%d/%m/%Y"))
fecpa$SamplingDate= format(as.Date(fecpa$SamplingDate, format="%d/%m/%Y"))
fecpa$DateOfBirth = format(as.Date(fecpa$DateOfBirth, format="%d/%m/%Y"))
fecpa$age = difftime(fecpa$AnalysisDate, fecpa$DateOfBirth, units = "days")

## Rename mare type
fecpa$Type = as.character(fecpa$Type)
fecpa$Type[fecpa$Type=='WS'] = 'MARE'
fecpa$Type = factor(fecpa$Type)

## Define families
family = fam[,c('MARE','FOAL')]
fecpa$fami = 0
fecpa$fami[fecpa$Type=='MARE'] = match(fecpa$Id[fecpa$Type=='MARE'],family$MARE)
fecpa$fami[fecpa$Type=='FOAL'] = match(fecpa$Id[fecpa$Type=='FOAL'],family$FOAL)
fecpa$fami = factor(fecpa$fami)

head(fecpa)
#      Id SamplingDate AnalysisDate eggsStr FECs eggsParascaris FECp  Date Pasture Type DateOfBirth
# 1 MW585   2019-06-18   2019-06-21       0    0              0    0 06/21      PA FOAL  2019-05-25
# 2 MW585   2019-07-17   2019-07-18       0    0              0    0 07/18      PA FOAL  2019-05-25
# 3 MW585   2019-08-19   2019-08-20       1   15              0    0 08/20      PA FOAL  2019-05-25
# 4 MW585   2019-09-16   2019-09-17       3   45              8  120 09/17      PA FOAL  2019-05-25
# 5 MW585   2019-08-30   2019-08-30      25  375              1   15 08/30      PA FOAL  2019-05-25
# 6 MW586   2019-07-17   2019-07-18       0    0              0    0 07/18      PA FOAL  2019-06-01
#        age fami
# 1  27 days    3
# 2  54 days    3
# 3  87 days    3
# 4 115 days    3
# 5  97 days    3
# 6  47 days    4

###Export
write.csv(fecpa,file = './supplementary_Table1.csv',quote=F,row.names=F)

###----------------------------============================== FEC trajectory ==================================-----------------------------

### Summary statistics for FEC in mares
aggregate(FECs ~ AnalysisDate, data = fecpa[fecpa$Type=='MARE',],FUN=summary)
# AnalysisDate  FECs.Min. FECs.1st Qu. FECs.Median  FECs.Mean FECs.3rd Qu.  FECs.Max.
# 1   2019-06-03    0.00000      0.00000     0.00000   34.61538      0.00000  400.00000
# 2   2019-07-18    0.00000      3.75000    30.00000  234.64286    341.25000 1410.00000
# 3   2019-08-20    0.00000     15.00000    75.00000  463.92857    386.25000 2520.00000
# 4   2019-08-30   30.00000     63.75000   202.50000  208.92857    296.25000  450.00000
# 5   2019-09-17    0.00000      0.00000     0.00000  141.42857    198.75000 1080.00000

aggregate(FECp ~ AnalysisDate, data = fecpa[fecpa$Type=='MARE',],FUN=summary)
# AnalysisDate FECp.Min. FECp.1st Qu. FECp.Median FECp.Mean FECp.3rd Qu. FECp.Max.
# 1   2019-06-03  0.000000     0.000000    0.000000  0.000000     0.000000  0.000000
# 2   2019-07-18  0.000000     0.000000    0.000000  0.000000     0.000000  0.000000
# 3   2019-08-20  0.000000     0.000000    0.000000  0.000000     0.000000  0.000000
# 4   2019-08-30  0.000000     0.000000    0.000000  2.142857     0.000000 30.000000
# 5   2019-09-17  0.000000     0.000000    0.000000  0.000000     0.000000  0.000000

### Age at 1st positive egg excretion in foals; strongyle vs. parascaris
aggregate(FECs ~ AnalysisDate, data = fecpa[fecpa$Type=='FOAL',],FUN=summary)
# AnalysisDate  FECs.Min. FECs.1st Qu. FECs.Median  FECs.Mean FECs.3rd Qu.  FECs.Max.
# 1   2019-06-21   0.00000      0.00000     0.00000  30.00000     37.50000 150.00000
# 2   2019-07-18   0.00000      0.00000     0.00000  11.78571      0.00000 120.00000
# 3   2019-08-20   0.00000      0.00000     7.50000  19.28571     15.00000  90.00000
# 4   2019-08-30   0.00000    105.00000   232.50000 227.14286    333.75000 570.00000
# 5   2019-09-17   0.00000     18.75000    37.50000  93.21429     45.00000 825.00000

aggregate(FECp ~ AnalysisDate, data = fecpa[fecpa$Type=='FOAL',],FUN=summary)
# AnalysisDate FECp.Min. FECp.1st Qu. FECp.Median FECp.Mean FECp.3rd Qu. FECp.Max.
# 1   2019-06-21   0.000000     0.000000    0.000000   0.000000     0.000000   0.000000
# 2   2019-07-18   0.000000     0.000000    0.000000   0.000000     0.000000   0.000000
# 3   2019-08-20   0.000000     0.000000    0.000000   0.000000     0.000000   0.000000
# 4   2019-08-30   0.000000     0.000000    0.000000   5.357143     0.000000  45.000000
# 5   2019-09-17   0.000000     0.000000    0.000000  32.142857     0.000000 210.000000

a = data.frame(aggregate(FECs ~ AnalysisDate + Type, data = fecpa, FUN=summary))
a
#    AnalysisDate Type FECs.Min. FECs.1st Qu. FECs.Median FECs.Mean FECs.3rd Qu. FECs.Max.
# 1    2019-06-21 FOAL      0.00         0.00        0.00     30.00        37.50    150.00
# 2    2019-07-18 FOAL      0.00         0.00        0.00     11.79         0.00    120.00
# 3    2019-08-20 FOAL      0.00         0.00        7.50     19.29        15.00     90.00
# 4    2019-08-30 FOAL      0.00       105.00      232.50    227.14       333.75    570.00
# 5    2019-09-17 FOAL      0.00        18.75       37.50     93.21        45.00    825.00
# 6    2019-06-03 MARE      0.00         0.00        0.00     34.62         0.00    400.00
# 7    2019-07-18 MARE      0.00         3.75       30.00    234.64       341.25   1410.00
# 8    2019-08-20 MARE      0.00        15.00       75.00    463.93       386.25   2520.00
# 9    2019-08-30 MARE     30.00        63.75      202.50    208.93       296.25    450.00
# 10   2019-09-17 MARE      0.00         0.00        0.00    141.43       198.75   1080.00

###----------------------------============================== IVM efficacy in mares ==================================-----------------------------
# W655 310KG 350 ERAQUELL dosis
# W662 316KG 350
# W675 265 KG 300
# W707 232 KG 250
# W688 424KG 450
# W650 318 KG 350
### Does IVM treatment in mares impact on foals egg excretion ?
drench.mare = c('W655', 'W662', 'W675', 'W707', 'W688', 'W650')
drench.foal = c('W762', 'MW591', 'MW590', 'MW586', 'W758','MW592')
drench = c(drench.mare,drench.foal)

## Define ivermectin datapoints
fecpa$Treatment = 'None'
fecpa$Treatment[fecpa$Id %in% drench & fecpa$AnalysisDate>'2019-08-29'] = 'Ivermectin'
fecpa$Treatment = factor(fecpa$Treatment)

df0 = fecpa[fecpa$AnalysisDate=='2019-08-30',c('Id','FECs','FECp','Pasture','Type','fami')]
colnames(df0)[2]='FEC0'
colnames(df0)[3]='FECp0'

df14 = fecpa[fecpa$AnalysisDate=='2019-09-17',c('Id','FECs','FECp')]
colnames(df14)[2]='FEC14'
colnames(df14)[3]='FECp14'

df = merge(df0,df14,by = 'Id')
df$GP = 'TM'
df$GP[df$Id %in% drench] = 'IVM'

head(df)
#      Id FEC0 FECp0 Pasture Type fami FEC14 FECp14  GP
# 1 MW585  375    15      PA FOAL    3    45    120  TM
# 2 MW586   30    15      PA FOAL    4    45      0 IVM
# 3 MW587  150     0      PA FOAL    1    30      0  TM
# 4 MW588  270     0      PA FOAL    9    45      0  TM
# 5 MW589  360     0      PA FOAL   14    15    120  TM
# 6 MW590   45     0      PN FOAL   10    30      0 IVM

aggregate(FEC0 ~ GP, FUN=mean,data=df) ## overall
#    GP   FEC0
# 1 IVM 223.75
# 2  TM 213.75

aggregate(FEC14 ~ GP, FUN=mean,data=df) ## overall
#    GP  FEC14
# 1 IVM  18.75
# 2  TM 191.25

aggregate(FEC0 ~ GP, FUN=mean,data=df[df$Type=='FOAL',]) ## Foals
#    GP FEC0
# 1 IVM  115
# 2  TM  311

aggregate(FEC14 ~ GP, FUN=mean,data=df[df$Type=='FOAL',]) ## Foals
#    GP FEC14
# 1 IVM  37.5
# 2  TM 135.0

### FECR without accounting for FEC trajectory in control individuals // overall

#- FEC14/FECJ0 post ivm
df$FECR[df$GP=='IVM'] = (1-df$FEC14[df$GP=='IVM']/df$FEC0[df$GP=='IVM'])*100
na.omit(df[df$GP=='IVM',])
#       Id FEC0 FECp0 Pasture Type FEC14 FECp14  GP      FECR
# 2  MW586   30    15      PA FOAL    45      0 IVM -50.00000
# 6  MW590   45     0      PN FOAL    30      0 IVM  33.33333
# 7  MW591  180     0      PN FOAL    45      0 IVM  75.00000
# 8  MW592   90    45      PA FOAL    75      0 IVM  16.66667
# 9   W650  435     0      PA MARE     0      0 IVM 100.00000
# 11  W655  300     0      PN MARE     0      0 IVM 100.00000
# 12  W662  360     0      PN MARE     0      0 IVM 100.00000
# 13  W675  450     0      PN MARE     0      0 IVM 100.00000
# 18  W688  165     0      PA MARE     0      0 IVM 100.00000
# 21  W707  285     0      PA MARE     0      0 IVM 100.00000
# 24  W758  345     0      PA FOAL    30      0 IVM  91.30435

summary(df$FECR[df$GP=='IVM'])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  -50.00   54.17  100.00   69.66  100.00  100.00       1 

summary(df$FECR[df$GP=='IVM' & df$Type=='FOAL'])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   -50.00   16.67   33.33   33.26   75.00   91.30       1 

summary(df$FECR[df$GP=='IVM' & df$Type=='MARE'])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100     100     100     100     100     100 

###--- FECR accounting for FEC trajectory in control group
#(FEC14/FECJ0)ivm x (FECJ0/FEC14)tm    
100*(1-mean(df$FEC14[df$GP=='IVM' & df$Type=='FOAL' & df$FEC0!=0]/df$FEC0[df$GP=='IVM' & df$Type=='FOAL' & df$FEC0!=0],na.rm=T)*
  mean(df$FEC0[df$GP=='TM' & df$Type=='FOAL' & df$FEC0!=0]/df$FEC14[df$GP=='TM' & df$Type=='FOAL' & df$FEC0!=0],na.rm=T))
#[1] -554

(1-mean(df$FEC14[df$GP=='IVM'& df$Type=='MARE']/df$FEC0[df$GP=='IVM' & df$Type=='MARE'],na.rm=T)*mean(df$FEC0[df$GP=='TM' & df$Type=='MARE']/df$FEC14[df$GP=='TM' & df$Type=='MARE'],na.rm=T))*100
#[1] 100

###--- FECR in foals considering the only individuals excreting more than 150 eggs/g
ivmf150 = df[df$GP=='IVM' & df$Type=='FOAL' & df$FEC0 >150,]
ivmf150
#       Id FEC0 FECp0 Pasture Type fami FEC14 FECp14  GP FECR
# 7  MW591  180     0      PN FOAL    6    45      0 IVM 75.0
# 24  W758  345     0      PA FOAL    8    30      0 IVM 91.3

# FECR
red = ivmf150$FEC14/ivmf150$FEC0
1-ivmf150$FEC14/ivmf150$FEC0
#[1] 0.750 0.913 

# With control group
tmf150 = df[df$GP=='TM' & df$Type=='FOAL' & df$FEC0 >150,]
tmf150
#       Id FEC0 FECp0 Pasture Type fami FEC14 FECp14 GP FECR
# 1  MW585  375    15      PA FOAL    3    45    120 TM   NA
# 4  MW588  270     0      PA FOAL    9    45      0 TM   NA
# 5  MW589  360     0      PA FOAL   14    15    120 TM   NA
# 23  W757  195     0      PN FOAL    5    15      0 TM   NA
# 25  W759  300     0      PN FOAL    2    90      0 TM   NA
# 26  W760  270     0      PN FOAL    7    15      0 TM   NA
# 27  W761  570     0      PN FOAL   11   825    210 TM   NA

corfac = tmf150$FEC0/tmf150$FEC14
corfac
#[1]  8.333  6.000 24.000 13.000  3.333 18.000  0.691
100*(1 - red*mean(corfac))
#[1] -161.99    8.87

## Confidence interval for FECR
library(eggCounts)
#options(digits=5)
b0.1 = fecrtCI(df$FEC0[df$GP=='IVM' & df$Type=='FOAL'],df$FEC14[df$GP=='IVM' & df$Type=='FOAL'],paired=TRUE,R=1000,alpha=.05)
b0.1
# $estimate
# [1] 67.3913
# 
# $bootCI
# [1]  0.00000 85.41455
# 
# $approxCI
# [1] -6.341105 90.000790

b0.2 = fecrtCI(df$FEC0[df$GP=='IVM' & df$Type=='FOAL' & df$FEC0>0],
               df$FEC14[df$GP=='IVM' & df$Type=='FOAL' & df$FEC0> 0],paired=TRUE,R=1000,alpha=.05)
b0.2
# $estimate
# [1] 67.3913
# 
# $bootCI
# [1]  6.25000 85.54217
# 
# $approxCI
# [1]  6.22162 88.66128

### Significant reduction in control group
b1.2 = fecrtCI(df$FEC0[df$GP=='TM' & df$Type=='FOAL' & df$FEC0>0],
               df$FEC14[df$GP=='TM' & df$Type=='FOAL' & df$FEC0> 0],paired=TRUE,R=1000,alpha=.05)
b1.2
# $estimate
# [1] 56.6
# 
# $bootCI
# [1] 10.4 90.8
# 
# $approxCI
# [1] -115.6   91.3

###--- Test whether observed difference is significant and relate to other variables (age, sex)
## Create table for the two considered time points (d0 and d14)
fecrdf = fecpa[fecpa$SamplingDate>='2019-08-30',]
fecrdf$Sex = 'F'
fecrdf$Sex = as.character(fecrdf$Sex)
fecrdf$Sex[fecrdf$Type=='FOAL'] = as.character(fam$SEX[match(fecrdf$Id[fecrdf$Type=='FOAL'],fam$FOAL)])
fecrdf$Sex = factor(fecrdf$Sex)
fecrdf$Treatment = factor(fecrdf$Treatment,levels = c('None','Ivermectin'))

## Nb individuals available within each treatment group
table(fecrdf$Treatment[fecrdf$Type=='FOAL'],fecrdf$Sex[fecrdf$Type=='FOAL'])/length(levels(as.factor(fecrdf$SamplingDate)))
#            F M
# None       4 4
# Ivermectin 2 4

qqcomp(fitdist(data = fecrdf$FECs, distr = "norm", method = "mle"),main="FEC \n QQ-Plot - Normal")
qqcomp(fitdist(data = fecrdf$FECs, distr = "pois", method = "mle"),main="FEC \n QQ-Plot - Poisson")
qqcomp(fitdist(data = fecrdf$FECs, distr = "nbinom", method = "mle"),main="FEC \n QQ-Plot - Negative binomial")

###-- Efficacy as a function of age: slightly better efficacy in older foals
dfoal = df[df$Type=='FOAL',]
dfoal$Age = fecrdf$age[match(dfoal$Id,fecrdf$Id)]

ggplot(dfoal,aes(x=Age, y = FECR))+
  geom_point() 

###--- GEE model
dfgee = dfoal[dfoal$Type=='FOAL',c('Id','GP','FEC0')]
colnames(dfgee)[3]='FEC'
dfgee$day = 0
tmp = dfoal[dfoal$Type=='FOAL',c('Id','GP','FEC14')]
colnames(tmp)[3]='FEC'
tmp$day = 1
dfgee = rbind(dfgee,tmp)
dfgee = dfgee[order(dfgee$Id,dfgee$day),]
dfgee$day = factor(dfgee$day)

#Remove W762 that has no egg
dfgee = dfgee[dfgee$Id!='W762',]
dfgee$Id = factor(dfgee$Id)
dfgee$GP = factor(dfgee$GP, levels=c('TM','IVM'))

a=aggregate(FEC ~ GP, data= dfgee, FUN=mean)
a$std = aggregate(FEC ~ GP, data= dfgee, FUN=sd)
a
#    GP   FEC std.GP std.FEC
# 1  TM 223.1     TM     229
# 2 IVM  91.5    IVM     100

a=aggregate(FEC ~ GP+day, data= dfgee, FUN=mean)
a$std = aggregate(FEC ~ GP+day, data= dfgee, FUN=sd)
a
#    GP day FEC std.GP std.day std.FEC
# 1  TM   0 311     TM       0   129.0
# 2 IVM   0 138    IVM       0   129.6
# 3  TM   1 135     TM       1   279.9
# 4 IVM   1  45    IVM       1    18.4

mf = formula(FEC ~ GP*day, data=dfgee) 

### Model FECR for each farm using neg binom
fit.full <- geem(FEC ~ GP*day,
                 family = negative.binomial(theta=1.265,link = 'log'), 
                 data = dfgee, 
                 id = Id, 
                 corstr = "ar1") ##correlation decreases as a power of how many timepoints apart 2 observations// eq to exch in this case as only 2 tp

summary(fit.full)
#             Estimates Model SE Robust SE   wald      p
# (Intercept)     4.927    0.568     0.376 13.110 0.0000
# GPTM            0.813    0.723     0.400  2.033 0.0420
# day1           -1.121    0.672     0.457 -2.454 0.0141
# GPTM:day1       0.285    0.855     0.737  0.387 0.6988
# 
# Estimated Correlation Parameter:  0.306 
# Correlation Structure:  ar1 
# Est. Scale Parameter:  2.02 
# 
# Number of GEE iterations: 2 
# Number of Clusters:  13    Maximum Cluster Size:  2 
# Number of observations with nonzero weight:  26

## Extract odds ratio
est = coef(fit.full)
s = summary(fit.full)
upr = est + 1.96*s$se.robust
lwr = est - 1.96*s$se.robust
OR.CI = exp(cbind(est,lwr,upr))
OR.CI
#                 est    lwr     upr
# (Intercept) 138.000 66.071 288.236
# GPTM          2.255  1.030   4.940
# day1          0.326  0.133   0.798
# GPTM:day1     1.330  0.314   5.643

OR.CI = data.frame(OR.CI)
OR.CI$variable = row.names(OR.CI)
OR.CI

ggplot(OR.CI,aes(x = variable, y = est)) +
  geom_point(size = 4) + 
  geom_errorbar(aes(ymin=lwr,ymax=upr,width=0)) +
  scale_y_log10() +
  geom_hline(yintercept=1,col="lightgrey") + 
  theme_classic() + coord_flip()

### Figure 1
### FEC trajectory in mares and foals, by family

## Plot
fecpa$drench = 'None'
fecpa$drench[fecpa$Id %in% drench]='Ivermectin'

pdf(file = './Supplementary_Figure1.pdf') #, width=10,height=10)
ggplot(fecpa, 
       aes(x = AnalysisDate, y = FECs + 1, col = Type, shape = Treatment, group = Id)) + 
  geom_hline(yintercept = 150, col='black', lty = 2) +
  geom_point(size = 4, alpha = .6) + geom_line() + 
  facet_wrap(~ fami, ncol = 3) +
  scale_color_manual(values = viridis_pal(option='D')(4)[c(1,3)]) +
  scale_y_log10(limits = c(1,3645), breaks = c(1, 15, 45, 135, 405, 1215,3645)) +
  ylab('Faecal Egg Count (eggs/g)') + xlab('Analysis date') +
  theme_classic() +
  theme(legend.position = 'bottom',
        text = element_text(size = 11),
        axis.text.x = element_text(angle=45, hjust=0.9, size =9))
dev.off()

dfgee$Timepoint = 'Before ivermectin treatment'
dfgee$Timepoint[dfgee$day==1] = 'Following ivermectin treatment (18 days)'
dfgee$Timepoint = factor(dfgee$Timepoint)
dfgee$Group = 'Control'
dfgee$Group[dfgee$GP=='IVM'] = 'Mare treated'

pdf(file = './Figure1.pdf') 
ggplot(dfgee, 
       aes(x = Timepoint, y = FEC+1, col = Group,group = Id)) +
  #geom_boxplot(alpha = .8) +
  geom_point(size = 3, alpha = .6) + geom_line() + 
  #scale_fill_manual(values = viridis_pal(option='D')(4)[c(1,3)]) +
  scale_color_manual(values = viridis_pal(option='D')(4)[c(1,3)]) +
  scale_y_log10() +#limits = c(1,1000), breaks = c(1, 3, 10, 30, 100, 300,1000)) +
  ylab('Faecal Egg Count (eggs/g)') + xlab('') +
  theme_classic() +
  theme(legend.position = 'bottom',
        text = element_text(size = 14),
        axis.text.x = element_text(angle=0, vjust=1, size = 12))
dev.off()

###----------------------------============================== IVM pharmacokinetics ==================================-----------------------------
ivmpk = read.csv(file = './ivm_pharmacokinetics.csv', sep = ';', header = T)

##get mares and foals (ordered by couples)
mares = sapply(stringr::str_split(colnames(ivmpk[2:7]),'_'),function(x) x[1])
foals = sapply(stringr::str_split(colnames(ivmpk[8:13]),'_'),function(x) x[1])

## refmt file
ivmpk = reshape2::melt(ivmpk)

##Define new variables
ivmpk$Hour = as.integer(as.character(gsub('T','',ivmpk$Time)))
ivmpk$Id = sapply(stringr::str_split(ivmpk$variable,'_'),function(x) x[1])
ivmpk$SampleType = sapply(stringr::str_split(ivmpk$variable,'_'),function(x) x[2])

## Retrieve pony status
type = unique(fecpa[,c('Id','Type')])
ivmpk = merge(ivmpk,type, by = 'Id')

head(ivmpk)
#      Id Time     variable value Hour SampleType Type
# 1 MW586   T1 MW586_Plasma    NA    1     Plasma FOAL
# 2 MW586   T8 MW586_Plasma  0.14    8     Plasma FOAL
# 3 MW586   T0 MW586_Plasma    NA    0     Plasma FOAL
# 4 MW586  T96 MW586_Plasma  0.30   96     Plasma FOAL
# 5 MW586 T168 MW586_Plasma  0.21  168     Plasma FOAL
# 6 MW586  T24 MW586_Plasma  0.42   24     Plasma FOAL

###--- Plot IVM pharmacokinetics
colist = c('#bf812d','#35978f','#8073ac')

## Plasma from mares
plm = ggplot(ivmpk[ivmpk$Type=='MARE' & ivmpk$SampleType=='Plasma',],aes(x = Hour, y = value, group = Id)) +
  geom_point(col = colist[1]) + geom_line(col = colist[1]) +
  ylab('Ivermectin concentration (ng/mL)') + xlab ('Hour after ivermectin treatment') +
  ggtitle('Mare plasma')

## Plasma from foals
plf = ggplot(ivmpk[ivmpk$Type=='FOAL' & ivmpk$SampleType=='Plasma',],aes(x = Hour, y = value, group = Id, col = colist[2])) +
  geom_point(col = colist[2]) + geom_line(col = colist[2]) +
  ylab('Ivermectin concentration (ng/mL)') + xlab ('Hour after ivermectin treatment') +
  ggtitle('Foal plasma')

## Milk from mares
mim = ggplot(ivmpk[ivmpk$Type=='MARE' & ivmpk$SampleType=='Milk',],aes(x = Hour, y = value, group = Id, col = colist[3])) +
  geom_point(col = colist[3]) + geom_line(col = colist[3]) +
  ylab('Ivermectin concentration (ng/mL)') + xlab ('Hour after ivermectin treatment') +
  ggtitle('Mare milk')

multiplot(plm,plf,mim)

## Couple plots between mare and their respective foals
couple = ivmpk[ivmpk$Type!='FOAL',]
couple$Id = factor(couple$Id)
dim(couple)
#[1] 84  7

ivmfoals = na.omit(ivmpk[ivmpk$Type=='FOAL',]) ## remove T0 and T+1h as they do not exist for foals
ivmfoals$Id = mares[match(ivmfoals$Id,foals)]
ivmfoals$SampleType = 'Foal plasma'

couple = rbind(couple,ivmfoals)
dim(couple)
#[1] 114   7

aggregate(value ~ Id,FUN=mean,data=ivmpk[ivmpk$Hour==8 & ivmpk$Type=='MARE',])
#     Id  value
# 1 W650 28.745
# 2 W655 26.640
# 3 W662 28.125
# 4 W675 19.635
# 5 W688 43.325
# 6 W707 18.685

pdf(file = './Figure2.pdf')
ggplot(couple,aes(x = Hour, y = value+1, group = paste0(Id,SampleType),shape = SampleType, col = SampleType)) +
  geom_point() + geom_line() + facet_wrap(~ Id) + theme_light() +
  scale_color_manual(values=cols_tw[c(3,9,6)]) + 
  scale_x_continuous(limits = c(0,170), breaks = seq(0, 170, 20)) +
  scale_y_log10(limits = c(1,81), breaks = c(0, 2, 8, 26, 80)) +
  ylab('Ivermectin concentration (ng/mL)') + xlab ('Hours post-treatment') +
  theme(legend.position = 'bottom',
        text = element_text(size = 12))
dev.off()

## Corr between mare and foal
milk = couple[grep('Milk',couple$variable),c('Id','Time','value')]
colnames(milk)[3]='milk'
mplasma = couple[couple$Type=='MARE' & couple$SampleType=='Plasma',c('Id','Time','value')]
colnames(mplasma)[3]='mplasma'
fplasma = couple[couple$Type=='FOAL' & couple$SampleType=='Foal plasma',c('Id','Time','value')]
colnames(fplasma)[3]='fplasma'

a=merge(milk, mplasma, by =c('Id', 'Time'))
b=merge(a,fplasma, by = c('Id','Time'))

rcorr(as.matrix(b[,3:5]), type ='spearman')
#          milk mplasma fplasma
# milk     1.00    0.61   -0.15
# mplasma  0.61    1.00   -0.12
# fplasma -0.15   -0.12    1.00
# 
# n= 30 
# 
# 
# P
#        milk   mplasma fplasma
# milk           0.0003  0.4340 
# mplasma 0.0003         0.5445 
# fplasma 0.4340 0.5445     

