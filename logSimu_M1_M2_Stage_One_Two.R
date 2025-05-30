rm(list=ls())
setwd('D:/project/sbr-template/TR_BBD_SBR - Copy/Suplementary_R_Code')

# main functions
source("homoBinomial.R")            # merging function
source("postProbBasketM1.R")        # method 1 two stage
source("postProbBasketM2.R")        # method 2 two stage
source("postProbBasketM1s1.R")      # method 1 single stage 
source("postProbBasketM2s1.R")      # method 2 single stage

#other functions to calibrate Ce and plot the PCS under different sample sizes
source("calibCEwithn1.R")           # generate Ce at fixed n1=c(10:30) and fixed n2=20 UNDER TWO-STAGE for both methods
source("calibCEwithn2.R")           # generate Ce at fixed n1=c(10:30) and fixed n2=20 UNDER TWO-STAGE for both methods
source("sampleSizPCSDiffn1.R")      # differrent n1=c(10:30) sample size vs PCS two stage at fixed n1=20
source("sampleSizPCSDiffn2.R")      # differrent n2=c(10:30) sample size vs PCS two stage at fixed n1=20
source("Ce.R")                      # for all cf, CeI and CeII from the given pi_{jk} 
source("mergeRate.R")               # used to check the merging rate of merging function homoBinomial.R


# pi_{jk}: prevalence of each cancer under each cutoff
PI=matrix(c(0.03,0.05,0.10,0.10,
            0.04,0.10,0.10,0.10,
            0.03,0.05,0.10,0.20),3,byrow=TRUE) 

#Here we display six scenario's based on different response rate p_{jk} (i.e. pM)

#Scenario 1
pM=rbind(c(0.20,0.20,0.45,0.45),
         c(0.20,0.20,0.45,0.45),
         c(0.20,0.20,0.45,0.45))
s1<-list( M1=postProbBasketM1(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.652, aP=1, bP=1),
          M2=postProbBasketM2(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.490,
                              row_cutEff=c(0.7716055, 0.7761992, 0.7617430),
                              col_cutEff=c(0.7742698, 0.7681103, 0.7755725, 0.7681103),
                              cutFisher=0.80, aP=1, bP=1),
          M1s1=postProbBasketM1s1(PI,pM,n=40,cutEff=0.761,
                                  aP=1, bP=1),
          M2s1=postProbBasketM2s1(PI,pM,n=40,cutEff=0.526, 
                                  row_cutEff=c(0.7984352, 0.8025440, 0.7889145),
                                  col_cutEff=c(0.7796536, 0.7759109, 0.7804498, 0.7759109),
                                  cutFisher=0.80, aP=1, bP=1))

#Scenario 2
pM=rbind(c(0.15,0.15,0.15,0.15),
         c(0.45,0.45,0.45,0.45),
         c(0.45,0.45,0.45,0.45))
s2<-list( M1=postProbBasketM1(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.652, aP=1, bP=1),
          M2=postProbBasketM2(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.490,
                              row_cutEff=c(0.7716055, 0.7761992, 0.7617430),
                              col_cutEff=c(0.7742698, 0.7681103, 0.7755725, 0.7681103),
                              cutFisher=0.80, aP=1, bP=1),
          M1s1=postProbBasketM1s1(PI,pM,n=40,cutEff=0.761,
                                  aP=1, bP=1),
          M2s1=postProbBasketM2s1(PI,pM,n=40,cutEff=0.526, 
                                  row_cutEff=c(0.7984352, 0.8025440, 0.7889145),
                                  col_cutEff=c(0.7796536, 0.7759109, 0.7804498, 0.7759109),
                                  cutFisher=0.80, aP=1, bP=1))

#Scenario 3
pM=rbind(c(0.20,0.20,0.20,0.45),
         c(0.20,0.20,0.45,0.45),
         c(0.20,0.45,0.45,0.45))
s3<-list( M1=postProbBasketM1(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.652, aP=1, bP=1),
          M2=postProbBasketM2(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.490,
                              row_cutEff=c(0.7716055, 0.7761992, 0.7617430),
                              col_cutEff=c(0.7742698, 0.7681103, 0.7755725, 0.7681103),
                              cutFisher=0.80, aP=1, bP=1),
          M1s1=postProbBasketM1s1(PI,pM,n=40,cutEff=0.761,
                                  aP=1, bP=1),
          M2s1=postProbBasketM2s1(PI,pM,n=40,cutEff=0.526, 
                                  row_cutEff=c(0.7984352, 0.8025440, 0.7889145),
                                  col_cutEff=c(0.7796536, 0.7759109, 0.7804498, 0.7759109),
                                  cutFisher=0.80, aP=1, bP=1))

#Scenario 4
pM=rbind(c(0.35,0.35,0.45,0.45),
         c(0.35,0.35,0.45,0.45),
         c(0.35,0.35,0.45,0.45))
s4<-list( M1=postProbBasketM1(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.652, aP=1, bP=1),
          M2=postProbBasketM2(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.490,
                              row_cutEff=c(0.7716055, 0.7761992, 0.7617430),
                              col_cutEff=c(0.7742698, 0.7681103, 0.7755725, 0.7681103),
                              cutFisher=0.80, aP=1, bP=1),
          M1s1=postProbBasketM1s1(PI,pM,n=40,cutEff=0.761,
                                  aP=1, bP=1),
          M2s1=postProbBasketM2s1(PI,pM,n=40,cutEff=0.526, 
                                  row_cutEff=c(0.7984352, 0.8025440, 0.7889145),
                                  col_cutEff=c(0.7796536, 0.7759109, 0.7804498, 0.7759109),
                                  cutFisher=0.80, aP=1, bP=1))

#Scenario 5
pM=rbind(c(0.35,0.35,0.35,0.45),
         c(0.35,0.35,0.45,0.45),
         c(0.35,0.45,0.45,0.45))
s5<-list( M1=postProbBasketM1(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.652, aP=1, bP=1),
          M2=postProbBasketM2(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.490,
                              row_cutEff=c(0.7716055, 0.7761992, 0.7617430),
                              col_cutEff=c(0.7742698, 0.7681103, 0.7755725, 0.7681103),
                              cutFisher=0.80, aP=1, bP=1),
          M1s1=postProbBasketM1s1(PI,pM,n=40,cutEff=0.761,
                                  aP=1, bP=1),
          M2s1=postProbBasketM2s1(PI,pM,n=40,cutEff=0.526, 
                                  row_cutEff=c(0.7984352, 0.8025440, 0.7889145),
                                  col_cutEff=c(0.7796536, 0.7759109, 0.7804498, 0.7759109),
                                  cutFisher=0.80, aP=1, bP=1))

#Scenario 6
pM=rbind(c(0.35,0.35,0.35,0.35),
         c(0.35,0.35,0.35,0.35),
         c(0.45,0.45,0.45,0.45))
s6<-list( M1=postProbBasketM1(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.652, aP=1, bP=1),
          M2=postProbBasketM2(PI,pM,n1=20,n2=20,alphaFut=0.20,cutEff=0.490,
                              row_cutEff=c(0.7716055, 0.7761992, 0.7617430),
                              col_cutEff=c(0.7742698, 0.7681103, 0.7755725, 0.7681103),
                              cutFisher=0.80, aP=1, bP=1),
          M1s1=postProbBasketM1s1(PI,pM,n=40,cutEff=0.761,
                                  aP=1, bP=1),
          M2s1=postProbBasketM2s1(PI,pM,n=40,cutEff=0.526, 
                                  row_cutEff=c(0.7984352, 0.8025440, 0.7889145),
                                  col_cutEff=c(0.7796536, 0.7759109, 0.7804498, 0.7759109),
                                  cutFisher=0.80, aP=1, bP=1))

# Output Summary 
PCS=list(Method1=c(s1$M1$PCS,s2$M1$PCS,s3$M1$PCS,s4$M1$PCS,s5$M1$PCS,s6$M1$PCS),
         Method2=c(s1$M2$PCS,s2$M2$PCS,s3$M2$PCS,s4$M2$PCS,s5$M2$PCS,s6$M2$PCS),
         Method1s1=c(s1$M1s1$PCS,s2$M1s1$PCS,s3$M1s1$PCS,s4$M1s1$PCS,s5$M1s1$PCS,s6$M1s1$PCS),
         Method2s1=c(s1$M2s1$PCS,s2$M2s1$PCS,s3$M2s1$PCS,s4$M2s1$PCS,s5$M2s1$PCS,s6$M2s1$PCS))

MSE=list(Method1=c(s1$M1$MSE,s2$M1$MSE,s3$M1$MSE,s4$M1$MSE,s5$M1$MSE,s6$M1$MSE),
         Method2=c(s1$M2$MSE,s2$M2$MSE,s3$M2$MSE,s4$M2$MSE,s5$M2$MSE,s6$M2$MSE),
         Method1s1=c(s1$M1s1$MSE,s2$M1s1$MSE,s3$M1s1$MSE,s4$M1s1$MSE,s5$M1s1$MSE,s6$M1s1$MSE),
         Method2s1=c(s1$M2s1$MSE,s2$M2s1$MSE,s3$M2s1$MSE,s4$M2s1$MSE,s5$M2s1$MSE,s6$M2s1$MSE))

PCSk=list(Method1=matrix(c(s1$M1$PCSk,s2$M1$PCSk,s3$M1$PCSk,s4$M1$PCSk,s5$M1$PCSk,s6$M1$PCSk),6,4, byrow = TRUE),
          Method2=matrix(c(s1$M2$PCSk,s2$M2$PCSk,s3$M2$PCSk,s4$M2$PCSk,s5$M2$PCSk,s6$M2$PCSk),6,4, byrow = TRUE),
          Method1s1=matrix(c(s1$M1s1$PCSk,s2$M1s1$PCSk,s3$M1s1$PCSk,s4$M1s1$PCSk,s5$M1s1$PCSk,s6$M1s1$PCSk),6,4, byrow = TRUE),
          Method2s1=matrix(c(s1$M2s1$PCSk,s2$M2s1$PCSk,s3$M2s1$PCSk,s4$M2s1$PCSk,s5$M2s1$PCSk,s6$M2s1$PCSk),6,4, byrow = TRUE))

R_saving=list(Method1=c(sum(s1$M1$R),sum(s2$M1$R),sum(s3$M1$R),sum(s4$M1$R),sum(s5$M1$R),sum(s6$M1$R)),
              Method2=c(sum(s1$M2$R),sum(s2$M2$R),sum(s3$M2$R),sum(s4$M2$R),sum(s5$M2$R),sum(s6$M2$R)))

efficacy_rate=list(Method1=list(s1$M1$PE,s2$M1$PE,s3$M1$PE,s4$M1$PE,s5$M1$PE,s6$M1$PE),
                   Method2=list(s1$M2$PE,s2$M2$PE,s3$M2$PE,s4$M2$PE,s5$M2$PE,s6$M2$PE),
                   Method1s1=list(s1$M1s1$PE,s2$M1s1$PE,s3$M1s1$PE,s4$M1s1$PE,s5$M1s1$PE,s6$M1s1$PE),
                   Method2s1=list(s1$M2s1$PE,s2$M2s1$PE,s3$M2s1$PE,s4$M2s1$PE,s5$M2s1$PE,s6$M2s1$PE))

nonFutile_rate=list(Method1=list(s1$M1$PnonFut,s2$M1$PnonFut,s3$M1$PnonFut,s4$M1$PnonFut,s5$M1$PnonFut,s6$M1$PnonFut),
                    Method2=list(s1$M2$PnonFut,s2$M2$PnonFut,s3$M2$PnonFut,s4$M2$PnonFut,s5$M2$PnonFut,s6$M2$PnonFut))
