# File: FLLargeRivers.R
# Purpose: Illustrate estimation for Florida Large lakes using SPlus
# Programmer: Tony Olsen
# re formatted for R
# edited by Chris Sedlacek
# on 15 August 2019

## Reanalysis done 10 24 2019 Jay Silvanima
## Exclusion file was incorrect and lat/longs were not formatted correctly in original file.


# Load psurvey.analysis library using menu
#  Packages: Load Packages  click on sp and spsurvey.  Also WriteXLS in order to export for NNc


# Read in data from Florida using File menu
LL.SITES<-LL_EXCLUSIONS_2016_2018_REDO_102419
names(LL.SITES)

# convert to Decimal degrees and do map projection
deg <- floor(LL.SITES$RANDOM_LATITUDE/10000)
min <- floor((LL.SITES$RANDOM_LATITUDE - deg*10000)/100)
sec <- LL.SITES$RANDOM_LATITUDE - deg*10000 - min*100
LL.SITES$latdd <- deg + min/60 + sec/3600
deg <- floor(LL.SITES$RANDOM_LONGITUDE/10000)
min <- floor((LL.SITES$RANDOM_LONGITUDE - deg*10000)/100)
sec <- LL.SITES$RANDOM_LONGITUDE - deg*10000 - min*100
LL.SITES$londd <- deg + min/60 + sec/3600

# do equal area projection for variance estimation
tmp <- marinus(LL.SITES$latdd, LL.SITES$londd)
LL.SITES$xmarinus <- tmp[,'x']
LL.SITES$ymarinus <- tmp[,'y']

# check design variables
table(LL.SITES$CAN_BE_SAMPLED)
table(LL.SITES$EXCLUSION_CATEGORY)
table(LL.SITES$EXCLUSION_CRITERIA)
table(LL.SITES$EXCLUSION_CATEGORY, LL.SITES$CAN_BE_SAMPLED)

# create status and TNT variables
LL.SITES$EXCLUSION_CATEGORY <- as.character(LL.SITES$EXCLUSION_CATEGORY)
LL.SITES$EXCLUSION_CATEGORY[LL.SITES$CAN_BE_SAMPLED == 'Y'] <- 'SAMPLED'
LL.SITES$EXCLUSION_CATEGORY <- as.factor(LL.SITES$EXCLUSION_CATEGORY)
levels(LL.SITES$EXCLUSION_CATEGORY)
LL.SITES$TNT <- LL.SITES$EXCLUSION_CATEGORY
levels(LL.SITES$TNT) <- list(T=c('SAMPLED', 'NO PERMISSION FROM OWNER', 'UNABLE TO ACCESS','OTHERWISE UNSAMPLEABLE','DRY'),
                             NT=c('WRONG RESOURCE/NOT PART OF TARGET POPULATION') )

table(LL.SITES$EXCLUSION_CATEGORY, LL.SITES$TNT)

##### adjust weights for use of over sample sites


tst <- !is.na(LL.SITES$TNT)
LL.SITES <- LL.SITES[tst,]

## Create basin variable
LL.SITES$basin <-  LL.SITES$REPORTING_UNIT
table(LL.SITES$basin)

framesize <- c("ZONE 1"=20352.681,"ZONE 2"=7659.159,
               "ZONE 3"=120941.938,"ZONE 4"=44684.968,"ZONE 5"=64552.710,
               "ZONE 6"=130039.903
) ### change to match stratum

nr <- nrow(LL.SITES)
LL.SITES$wgt <- adjwgt(rep(TRUE,nr), 
                       LL.SITES$NEST1_WT, 
                       LL.SITES$basin, 
                       framesize=framesize               
)							
tapply(LL.SITES$wgt, LL.SITES$basin, sum)


sites <- data.frame(siteID=LL.SITES$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=LL.SITES$PK_RANDOM_SAMPLE_LOCATION,
                     Combined=rep("All Basins", nr), basin=LL.SITES$basin)  ##rep('LL',nr) 
dsgn <- data.frame(siteID=LL.SITES$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=LL.SITES$wgt,
                   xcoord=LL.SITES$xmarinus ,
                   ycoord=LL.SITES$ymarinus ,
                   stratum=LL.SITES$basin ##rep('LL',nr)
)
data.cat <- data.frame(siteID=LL.SITES$PK_RANDOM_SAMPLE_LOCATION, 
                       EXCLUSION.CATEGORY=LL.SITES$EXCLUSION_CATEGORY,
                       TNTStatus=LL.SITES$TNT
)


software.program <- 'R Studio'
ExtentEst <- cat.analysis(sites, subpop, dsgn, data.cat,
                          popsize=list(Combined=list("All Basins"=framesize),
                                       basin=list("ZONE 1"=framesize,"ZONE 2"=framesize,
                                                  "ZONE 3"=framesize,"ZONE 4"=framesize,"ZONE 5"=framesize,
                                                  "ZONE 6"=framesize		
                                       )),
                          conf=95
)

### write out or export results
write.csv(ExtentEst,file='ExtentEst.csv')


###### Do summaries for data
# import data
Combined.LL<-LL_RESULTS_2016_2018

names(Combined.LL)
# look at sample type
summary(Combined.LL$Sample.Type)

# Note that have BLANK, BOTTOM and PRIMARY data in original spreadsheet.
# hence need to subset before do analyses
keep <- Combined.LL$Sample.Type == 'PRIMARY'

# merge with status
comb <- merge(LL.SITES, Combined.LL[keep,], by.x='PK_RANDOM_SAMPLE_LOCATION',
              by.y='Station.Name')

# check that have only PRIMARY data
summary(comb$Sample.Type)  # yes!
summary(comb$Matrix)
keep1 <- comb$Matrix == 'WATER'
comb<-(comb[keep1,])
summary(comb$Matrix)


########UIA Calculation####################

## Calculates unionized ammonia from total [NH4+NH3]

## Does not make ionic strength correction yet...

## Script requires temperature in degC, pH, and total ammonia value, in mg/L
## Warning! be sure to check odd UIA values, to see that required input variables exist
## Also, values are taken AS reported.  If total ammonia/ammonium is below detection, so will be the UIA!!!

## get temperature vector, convert to degK
vTemp<-comb$Water.Temperature
vTempK<-vTemp+273.2

## calculate equilibrium constant at sample Temperature
vEQconstant<-2729.92/vTempK + 0.0901821

## get UI fraction of total
vExponent<-vEQconstant - comb$pH..Field
vDenom<-10^vExponent+1
vUIAfraction<-1/vDenom

## multiply fraction times total for UIA in lab units
vUIAraw<-vUIAfraction*comb$Ammonia..Total..as.N

## now convert to criterion units and write to file
comb$UIA<-vUIAraw*1.2158

#########End UIA Calculation#################

######################################################################################
###################2006 305b TSI calculations#########################################
######################################################################################


## calculates tsi (1996)

## the following lines calculate vectors for the tsi variables and add them to the 

## data set as columns at the end of the file

comb$tntp<-(comb$Nitrate.Nitrite..Total..as.N+comb$Kjeldahl.Nitrogen..Total..as.N)/comb$Phosphorus..Total..as.P

comb$tsichl<-16.8+(14.4*log(comb$Chlorophyll.A..Monochromatic))

comb$tsitn<-56+(19.8*log(comb$Nitrate.Nitrite..Total..as.N+comb$Kjeldahl.Nitrogen..Total..as.N))

comb$tsitn2<-10*(5.96+2.15*log(comb$Nitrate.Nitrite..Total..as.N+comb$Kjeldahl.Nitrogen..Total..as.N+0.0001))

comb$tsitp<-(18.6*log(comb$Phosphorus..Total..as.P*1000))-18.4

comb$tsitp2<-10*(2.36*log(comb$Phosphorus..Total..as.P*1000)-2.38) 



## calculates tsi-nutr, for neither N nor P limiting

comb$tsinutr<-(comb$tsitp+comb$tsitn)/2



## calculate tsi-nutr for (1) P-limiting condition and (2) N-limiting condition

comb$tsinutr[comb$tntp > 30]<-comb$tsitp2[comb$tntp > 30]

comb$tsinutr[comb$tntp < 10]<-comb$tsitn2[comb$tntp < 10]   



## now average tsinutr and tsichl, doesn't that make sense? 

comb$XTSI<-(comb$tsichl+comb$tsinutr)/2



## now add another column, this time a tsi category, 1 for good, 2 for poor

## category is color-dependent, for low color 

comb$tsicat[comb$Color..true <= 40]<-1

comb$tsicat[comb$Color..true <= 40 & comb$XTSI > 40] <-2



## and for high color

comb$tsicat[comb$Color..true > 40]<-1

comb$tsicat[comb$Color..true > 40 & comb$XTSI > 60] <-2


##### Example continuous indicator estimation
nr <- nrow(comb)
sites <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins", nr), Basin=comb$basin )
dsgn <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb$wgt,
                   xcoord=comb$xmarinus ,
                   ycoord=comb$ymarinus ,
                   stratum=comb$basin ##rep('Combined.LL.04.08',nr)
)

data.cont<-data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION,comb[, seq(46, 68, by=2)] )
data.cont2 <-data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION,comb[, seq(72,122 , by=2)] )

ExampleEst <- cont.analysis(sites, subpop, dsgn, data.cont,
                            popsize=list(CombinedBasins=list("All Basins"=framesize),
                                         Basin=list("ZONE 1"=framesize,"ZONE 2"=framesize,
                                                    "ZONE 3"=framesize,"ZONE 4"=framesize,"ZONE 5"=framesize,
                                                    "ZONE 6"=framesize      
                                         ))
)
# write out the results
write.csv(ExampleEst$CDF,file='ExampleEstCDF.csv')
write.csv(ExampleEst$Pct,file='ExampleEstPCT.csv')
write.csv(ExampleEst$Tot,file='ExampleEstTOT.csv')

#ExampleEst <- cont.analysis(sites, subpop, dsgn, data.cont)
ExampleEst2 <- cont.analysis(sites, subpop, dsgn, data.cont2,
                             popsize=list(CombinedBasins=list("All Basins"=framesize),
                                          Basin=list("ZONE 1"=framesize,"ZONE 2"=framesize,
                                                     "ZONE 3"=framesize,"ZONE 4"=framesize,"ZONE 5"=framesize,
                                                     "ZONE 6"=framesize            
                                          ))
)




# write out the results for 2
write.csv(ExampleEst2$CDF,file='ExampleEstCDF2.csv')
write.csv(ExampleEst2$Pct,file='ExampleEstPCT2.csv')
write.csv(ExampleEst2$Tot,file='ExampleEstTOT2.csv')


######################################################################################
###################2006 305b Pie Charts (Lakes)#######################################
######################################################################################

###___________________________________________________________________________________

### Script to do categories for pie chartsUIA
UIAcat <- cut(comb$UIA, breaks=c(0,0.02,1000000), include.lowest=TRUE)
comb$UIAcat <- UIAcat
comb$UIAcat <- as.factor(comb$UIAcat)


nr <- nrow(comb)
sites <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb$basin) 
dsgn <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb$wgt,
                   xcoord=comb$xmarinus ,
                   ycoord=comb$ymarinus ,
                   stratum=comb$basin
)
data.catUIA <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                          UIAcategory=comb$UIAcat
)

CategoryEstUIA <- cat.analysis(sites, subpop, dsgn, data.catUIA,
                               popsize=list(CombinedBasins=list("All Basins"=framesize),
                                            Basin=list("ZONE 1"=framesize,"ZONE 2"=framesize,
                                                       "ZONE 3"=framesize,"ZONE 4"=framesize,"ZONE 5"=framesize,
                                                       "ZONE 6"=framesize          
                                            ))
)

# write out the results
write.csv(CategoryEstUIA,file='CategoryExampleUIA.csv')
###___________________________________________________________________________________

### Script to do categories for pie chartstsicat
tsicatcat <- cut(comb$tsicat, breaks=c(0,1,100), include.lowest=TRUE)
comb$tsicatcat <- tsicatcat
comb$tsicatcat <- as.factor(comb$tsicatcat)


nr <- nrow(comb)
sites <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb$basin) 
dsgn <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb$wgt,
                   xcoord=comb$xmarinus ,
                   ycoord=comb$ymarinus ,
                   stratum=comb$basin
)
data.cattsicat <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                             tsicatcategory=comb$tsicatcat
)

CategoryEsttsicat <- cat.analysis(sites, subpop, dsgn, data.cattsicat,
                                  popsize=list(CombinedBasins=list("All Basins"=framesize),
                                               Basin=list("ZONE 1"=framesize,"ZONE 2"=framesize,
                                                          "ZONE 3"=framesize,"ZONE 4"=framesize,"ZONE 5"=framesize,
                                                          "ZONE 6"=framesize                   ))
)

# write out the results
write.csv(CategoryEsttsicat,file='CategoryExampletsicatV2.csv')

####
###___________________________________________________________________________________
### Script to do categories for pie chartsE.coli

## comb<-comb[sapply(comb[,"Escherichia.Coli.Quanti.Tray_VQ"],function(x)any(x==c('','NA','A','I','T','K','N','L','J','W','U','Z','Q','ABQ','AKQ','AQ','B','BQ','KQ','UQ'))),]

X31616cat <- cut(comb$Escherichia.Coli.Quanti.Tray, breaks=c(0,409,10000000), include.lowest=TRUE)
comb$X31616cat <- X31616cat
comb$X31616cat <- as.factor(comb$X31616cat)


nr <- nrow(comb)
sites <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb$basin) 
dsgn <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb$wgt,
                   xcoord=comb$xmarinus ,
                   ycoord=comb$ymarinus ,
                   stratum=comb$basin
)
data.cat31616 <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                            X31616category=comb$X31616cat
)

CategoryEst31616 <- cat.analysis(sites, subpop, dsgn, data.cat31616,
                                 popsize=list(CombinedBasins=list("All Basins"=framesize),
                                              Basin=list("ZONE 1"=framesize,"ZONE 2"=framesize,
                                                         "ZONE 3"=framesize,"ZONE 4"=framesize,"ZONE 5"=framesize,
                                                         "ZONE 6"=framesize             
                                              ))
)

# write out the results
write.csv(CategoryEst31616,file='CategoryExampleEcoli.csv')


###___________________________________________________________________________________

###___________________________________________________________________________________

### Script to do categories for pie charts406

comb <- merge(LL.SITES, Combined.LL.04.08[keep,], by.x='PK.RANDOM.SAMPLE.LOCATION',
              by.y='Station.Name')
### comb<-comb[sapply(comb[,"pH..Field.VQ"],function(x)any(x==c('','NA','A','I','T','K','N','L','J','W','U','Z'))),]

X406cat <- cut(comb$pH..Field, breaks=c(0,5.999,8.5,14), include.lowest=TRUE)
comb$X406cat <- X406cat
comb$X406cat <- as.factor(comb$X406cat)


nr <- nrow(comb)
sites <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb$basin) 
dsgn <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb$wgt,
                   xcoord=comb$xmarinus ,
                   ycoord=comb$ymarinus ,
                   stratum=comb$basin
)
data.cat406 <- data.frame(siteID=comb$PK_RANDOM_SAMPLE_LOCATION, 
                          X406category=comb$X406cat
)

CategoryEst406 <- cat.analysis(sites, subpop, dsgn, data.cat406,
                               popsize=list(CombinedBasins=list("All Basins"=framesize),
                                            Basin=list("ZONE 1"=framesize,"ZONE 2"=framesize,
                                                       "ZONE 3"=framesize,"ZONE 4"=framesize,"ZONE 5"=framesize,
                                                       "ZONE 6"=framesize           
                                            ))
)

# write out the results
write.csv(CategoryEst406,file='CategoryExample406pH.csv')

###___________________________________________________________________________________
###___________________________________________________________________________________
# need to add the min and max NNC values

comb3<-comb

comb3$Colcat<- ifelse((comb3$Color..true > 40) ,0,1)
comb3$Alkcat<- ifelse((comb3$Alkalinity..Total..as.CaCO3 > 20) ,0,1)

comb3$col.alk.cat<-paste(comb3$Colcat, comb3$Alkcat)

## Use new col.alk.cat character variable to assign thresholds for TN and TP
comb3$TN.Min<- ifelse(comb3$col.alk.cat=="0 0",1.27,
                      ifelse(comb3$col.alk.cat=="0 1", 1.27,
                             ifelse(comb3$col.alk.cat== "1 0", 1.05,
                                    ifelse(comb3$col.alk.cat=="1 1", 0.51,NA))))
comb3$TN.Max<- ifelse(comb3$col.alk.cat=="0 0",2.23,
                      ifelse(comb3$col.alk.cat=="0 1", 2.23,
                             ifelse(comb3$col.alk.cat== "1 0", 1.91,
                                    ifelse(comb3$col.alk.cat=="1 1", 0.93,NA))))
comb3$TP.Min<- ifelse((comb3$Colcat==0 & comb3$BRREG_NAME=="WEST CENTRAL"),0.49,
                      ifelse(comb3$col.alk.cat=="0 0",0.05,
                      ifelse(comb3$col.alk.cat=="0 1",0.05,
                             ifelse(comb3$col.alk.cat=="1 0", 0.03,
                                    ifelse(comb3$col.alk.cat=="1 1",0.01,NA)))))
comb3$TP.Max<- ifelse((comb3$Colcat==0 & comb3$BRREG_NAME=="WEST CENTRAL"),0.49,
                      ifelse(comb3$col.alk.cat=="0 0",0.16,
                             ifelse(comb3$col.alk.cat=="0 1",0.16,
                                    ifelse(comb3$col.alk.cat=="1 0", 0.09,
                                           ifelse(comb3$col.alk.cat=="1 1",0.03,NA)))))



names(comb3)

comb3$TN<-(comb3$Kjeldahl.Nitrogen..Total..as.N+comb3$Nitrate.Nitrite..Total..as.N)

### Pass=1 AND Fail=0
## Using the minimum values

comb3$TNcat<-ifelse((comb3$TN.Min >= comb3$TN),1,0) 

comb3$TPcat<-ifelse((comb3$TP.Min >= comb3$Phosphorus..Total..as.P),1,0)  	
comb3$DOcat<-ifelse((comb3$DO.Conc >= comb3$Oxygen..Dissolved.Percent.Saturation),0,1) 

comb3$Tot<-(comb3$TNcat+comb3$TPcat+comb3$DOcat)


######################################################################################
###################2006 305b Pie Charts (Large Lakes)#######################################
######################################################################################


###___________________________________________________________________________________
### Script to do pie charts on overall total water quality
comb3$Totcat<-comb3$Tot
Totcat <- cut(comb3$Totcat, breaks=c(0,0.99,1.99,2.99,5), include.lowest=TRUE)
comb3$Totcat <- Totcat
comb3$Totcat <- as.factor(comb3$Totcat)


nr <- nrow(comb3)
sites <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb3$basin) 
dsgn <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb3$wgt,
                   xcoord=comb3$xmarinus ,
                   ycoord=comb3$ymarinus ,
                   stratum=comb3$basin
)
data.catTot <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                          Totcategory=comb3$Totcat
)

CategoryEstTot <-  cat.analysis(sites, subpop, dsgn, data.catTot,
                                popsize=list(CombinedBasins=list("All Basins"=framesize),
                                             Basin=list("ZONE 1"=framesize,
                                                        "ZONE 2"=framesize,
                                                        "ZONE 3"=framesize, 
                                                        "ZONE 4"=framesize,
                                                        "ZONE 5"=framesize,
                                                        "ZONE 6"=framesize 		
                                             ))
)
# write out the results
write.csv(CategoryEstTot,file='CategoryExampleToticatmin_LL.csv')

###
###_______________________________________________________________________________________


###___________________________________________________________________________________
### Script to do categories for pie charts of Total Nitrogen
comb3$TNicat<-comb3$TNcat
TNicat <- cut(comb3$TNicat, breaks=c(0,0.99,1), include.lowest=TRUE)
comb3$TNicat <- TNicat
comb3$TNicat <- as.factor(comb3$TNicat)


nr <- nrow(comb3)
sites <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb3$basin) 
dsgn <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb3$wgt,
                   xcoord=comb3$xmarinus ,
                   ycoord=comb3$ymarinus ,
                   stratum=comb3$basin
)
data.catTN <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                         TNcategory=comb3$TNicat
)

CategoryEstTN <-  cat.analysis(sites, subpop, dsgn, data.catTN,
                               popsize=list(CombinedBasins=list("All Basins"=framesize),
                                            Basin=list("ZONE 1"=framesize,
                                                       "ZONE 2"=framesize,
                                                       "ZONE 3"=framesize, 
                                                       "ZONE 4"=framesize,
                                                       "ZONE 5"=framesize,
                                                       "ZONE 6"=framesize		
                                            ))
)
# write out the results
write.csv(CategoryEstTN,file='CategoryExampleTNmin_LL.csv')

###
###_______________________________________________________________________________________
###
### Script to do categories for pie charts
comb3$TPicat<-comb3$TPcat
TPicat <- cut(comb3$TPicat, breaks=c(0,0.99,1), include.lowest=TRUE)
comb3$TPicat <- TPicat
comb3$TPicat <- as.factor(comb3$TPicat)


nr <- nrow(comb3)
sites <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb3$basin) 
dsgn <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb3$wgt,
                   xcoord=comb3$xmarinus ,
                   ycoord=comb3$ymarinus ,
                   stratum=comb3$basin
)

data.catTP <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                         TPcategory=comb3$TPicat
)

CategoryEstTP <-  cat.analysis(sites, subpop, dsgn, data.catTP,
                               popsize=list(CombinedBasins=list("All Basins"=framesize),
                                            Basin=list("ZONE 1"=framesize,
                                                       "ZONE 2"=framesize,
                                                       "ZONE 3"=framesize, 
                                                       "ZONE 4"=framesize,
                                                       "ZONE 5"=framesize,
                                                       "ZONE 6"=framesize 		
                                            ))
)
# write out the results
write.csv(CategoryEstTP,file='CategoryExampleTPmin_LL.csv')

###_______________________________________________________________________________________
###
### Script to do categories for pie charts
comb3$DOicat<-comb3$DOcat
DOicat <- cut(comb3$DOicat, breaks=c(0,0.99,1), include.lowest=TRUE)
comb3$DOicat <- DOicat
comb3$DOicat <- as.factor(comb3$DOicat)


nr <- nrow(comb3)
sites <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb3$basin) 
dsgn <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb3$wgt,
                   xcoord=comb3$xmarinus ,
                   ycoord=comb3$ymarinus ,
                   stratum=comb3$basin
)

data.catDO <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                         TPcategory=comb3$DOicat
)

CategoryEstDO <-  cat.analysis(sites, subpop, dsgn, data.catDO,
                               popsize=list(CombinedBasins=list("All Basins"=framesize),
                                            Basin=list("ZONE 1"=framesize,
                                                       "ZONE 2"=framesize,
                                                       "ZONE 3"=framesize, 
                                                       "ZONE 4"=framesize,
                                                       "ZONE 5"=framesize,
                                                       "ZONE 6"=framesize 		
                                            ))
)
# write out the results
write.csv(CategoryEstDO,file='CategoryExampleDOicat_LL.csv')

###_________________________________________________________________________________________

## Using the maximum values
comb3$TN2cat<-ifelse((comb3$TN.Max > comb3$TN),1,0) 

comb3$TP2cat<-ifelse((comb3$TP.Max >= comb3$Phosphorus..Total..as.P),1,0)  # comb2$TPicat<-0 else if (comb2$NNC.TP <comb2$Phosphorus..Total..as.P.) comb2$TPicat<-1

comb3$DOcat<-ifelse((comb3$DO.Conc > comb3$Oxygen..Dissolved.Percent.Saturation),0,1) 

comb3$Tot2<-(comb3$TN2cat+comb3$TP2cat+comb3$DOcat)


######################################################################################
###################2006 305b Pie Charts (LAKES)#######################################
######################################################################################


###___________________________________________________________________________________
### Script to do pie charts on overall total water quality
comb3$Tot2cat<-comb3$Tot2
Tot2cat <- cut(comb3$Tot2cat, breaks=c(0,0.99,1.99,2.99,5), include.lowest=TRUE)
comb3$Tot2cat <- Tot2cat
comb3$Tot2cat <- as.factor(comb3$Tot2cat)


nr <- nrow(comb3)
sites <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb3$basin) 
dsgn <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb3$wgt,
                   xcoord=comb3$xmarinus ,
                   ycoord=comb3$ymarinus ,
                   stratum=comb3$basin
)
data.catTot2 <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                           Totcategory=comb3$Tot2cat
)

CategoryEstTot2 <-  cat.analysis(sites, subpop, dsgn, data.catTot2,
                                 popsize=list(CombinedBasins=list("All Basins"=framesize),
                                              Basin=list("ZONE 1"=framesize,
                                                         "ZONE 2"=framesize,
                                                         "ZONE 3"=framesize, 
                                                         "ZONE 4"=framesize,
                                                         "ZONE 5"=framesize,
                                                         "ZONE 6"=framesize 		
                                              ))
)
# write out the results
write.csv(CategoryEstTot2,file='CategoryExampleTot2icatmax.csv')

###
###_______________________________________________________________________________________


###___________________________________________________________________________________
### Script to do categories for pie charts of Total Nitrogen
comb3$TN2icat<-comb3$TN2cat
TN2icat <- cut(comb3$TN2icat, breaks=c(0,0.99,1), include.lowest=TRUE)
comb3$TN2icat <- TN2icat
comb3$TN2icat <- as.factor(comb3$TN2icat)


nr <- nrow(comb3)
sites <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb3$basin) 
dsgn <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb3$wgt,
                   xcoord=comb3$xmarinus ,
                   ycoord=comb3$ymarinus ,
                   stratum=comb3$basin
)
data.catTN2 <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                          TNcategory=comb3$TN2icat
)

CategoryEstTN2 <-  cat.analysis(sites, subpop, dsgn, data.catTN2,
                                popsize=list(CombinedBasins=list("All Basins"=framesize),
                                             Basin=list("ZONE 1"=framesize,
                                                        "ZONE 2"=framesize,
                                                        "ZONE 3"=framesize, 
                                                        "ZONE 4"=framesize,
                                                        "ZONE 5"=framesize,
                                                        "ZONE 6"=framesize 		
                                             ))
)
# write out the results
write.csv(CategoryEstTN2,file='CategoryExampleTNmax_LL.csv')

###
###_______________________________________________________________________________________
###
### Script to do categories for pie charts
comb3$TP2icat<-comb3$TP2cat
TP2icat <- cut(comb3$TP2icat, breaks=c(0,0.99,1), include.lowest=TRUE)
comb3$TP2icat <- TP2icat
comb3$TP2icat <- as.factor(comb3$TP2icat)


nr <- nrow(comb3)
sites <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb3$basin) 
dsgn <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb3$wgt,
                   xcoord=comb3$xmarinus ,
                   ycoord=comb3$ymarinus ,
                   stratum=comb3$basin
)

data.catTP2 <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                          TPcategory=comb3$TP2icat
)

CategoryEstTP2 <-  cat.analysis(sites, subpop, dsgn, data.catTP2,
                                popsize=list(CombinedBasins=list("All Basins"=framesize),
                                             Basin=list("ZONE 1"=framesize,
                                                        "ZONE 2"=framesize,
                                                        "ZONE 3"=framesize, 
                                                        "ZONE 4"=framesize,
                                                        "ZONE 5"=framesize,
                                                        "ZONE 6"=framesize 		
                                             ))
)
# write out the results
write.csv(CategoryEstTP2,file='CategoryExampleTPmax_LL.csv')

###_______________________________________________________________________________________
###
### Script to do categories for pie charts
comb3$DO2icat<-comb3$DOcat
DO2icat <- cut(comb3$DO2icat, breaks=c(0,0.99,1), include.lowest=TRUE)
comb3$DO2icat <- DO2icat
comb3$DO2icat <- as.factor(comb3$DO2icat)


nr <- nrow(comb3)
sites <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb3$basin) 
dsgn <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb3$wgt,
                   xcoord=comb3$xmarinus ,
                   ycoord=comb3$ymarinus ,
                   stratum=comb3$basin
)

data.catDO2 <- data.frame(siteID=comb3$PK_RANDOM_SAMPLE_LOCATION, 
                          TPcategory=comb3$DO2icat
)

CategoryEstDO2 <-  cat.analysis(sites, subpop, dsgn, data.catDO2,
                                popsize=list(CombinedBasins=list("All Basins"=framesize),
                                             Basin=list("ZONE 1"=framesize,
                                                        "ZONE 2"=framesize,
                                                        "ZONE 3"=framesize, 
                                                        "ZONE 4"=framesize,
                                                        "ZONE 5"=framesize,
                                                        "ZONE 6"=framesize 		
                                             ))
)
# write out the results
write.csv(CategoryEstDO2,file='CategoryExampleDO2icat_LL.csv')
###
###_______________________________________________________________________________________
###
comb4<-comb3
names(comb4)

comb4$Chl.conc<- ifelse(comb4$col.alk.cat=="0 0",20,
                        ifelse(comb4$col.alk.cat=="0 1", 20,
                               ifelse(comb4$col.alk.cat== "1 0", 20,
                                      ifelse(comb4$col.alk.cat=="1 1", 6,NA))))

comb4$Chlpf<-ifelse((comb4$Chl.conc > comb4$Chlorophyll.A..Monochromatic),1,0) 

X32211cat <- cut(comb4$Chlpf, breaks=c(0,0.99,1), include.lowest=TRUE)
comb4$X32211cat <- X32211cat
comb4$X32211cat <- as.factor(comb4$X32211cat)


nr <- nrow(comb4)
sites <- data.frame(siteID=comb4$PK_RANDOM_SAMPLE_LOCATION, Use=rep(TRUE,nr) )
subpop <- data.frame(siteID=comb4$PK_RANDOM_SAMPLE_LOCATION, 
                     CombinedBasins=rep("All Basins",nr),
                     Basin=comb4$basin) 
dsgn <- data.frame(siteID=comb4$PK_RANDOM_SAMPLE_LOCATION, 
                   wgt=comb4$wgt,
                   xcoord=comb4$xmarinus ,
                   ycoord=comb4$ymarinus ,
                   stratum=comb4$basin
)
data.cat32211 <- data.frame(siteID=comb4$PK_RANDOM_SAMPLE_LOCATION, 
                            X32211category=comb4$X32211cat
)

CategoryEst32211 <- cat.analysis(sites, subpop, dsgn, data.cat32211,
                                 popsize=list(CombinedBasins=list("All Basins"=framesize),
                                              Basin=list("ZONE 1"=framesize ,"ZONE 2"=framesize ,
                                                         "ZONE 3"=framesize ,"ZONE 4"=framesize ,
                                                         "ZONE 5"=framesize ,"ZONE 6"=framesize	)
                                 ) )

# write out the results
write.csv(CategoryEst32211,file='Chlorophylla_pass_fail.csv')

###___________________________________________________________________________________