calibCEwithn1<-function(PI,n2=20,q0=0.20, q1=0.35,alphaFut=0.20,alphaEff=0.30,cutFisher=0.80,numSimu=1000000) 
{ # generate Ce at fixed n1=c(10:30) and fixed n2=20 UNDER TWO-STAGE
  n1=seq(10,30,1)
  JK=dim(PI); J=JK[1]; K=JK[2]
  Cf_n1=rep(0,length(n1)) 
  cutEff_In1=array(0,dim=c(length(n1),J,K));  cutEff_IIn1=array(0,dim=c(length(n1),J,K))                        
  row_cutEff=matrix(0,length(n1),J);   col_cutEff=matrix(0,length(n1),K)
  for (i in 1:length(n1)) { # length(n1)
    
    Ce_n1i=CeIIprof(n1=n1[i],n2,q0, q1,alphaFut,alphaEff,cutFisher,numSimu) #  all methods calibrated Cf and Ce for each n1 and fixed n2
    # threshold outputs
    Cf_n1[i]=Ce_n1i$cutFut          # calibrated Cf for each n1 at fixed n2
    cutEff_In1[i,,]=Ce_n1i$CeI      # calibrated CeI for each n1 at fixed n2
    row_cutEff[i,]=Ce_n1i$rowCeII   # calibrated rowCeII for each n1 at fixed n2
    col_cutEff[i,]=Ce_n1i$colCeII   # calibrated colCeII for each n1 at fixed n2
    cutEff_IIn1[i,,]=Ce_n1i$CeII    # calibrated CeII for each n1 at fixed n2
  }
  out<-list(Cf_n1=Cf_n1,cutEff_In1=cutEff_In1,
            row_cutEff=row_cutEff,col_cutEff=col_cutEff,
            cutEff_IIn1=cutEff_IIn1,
            n1=n1,n2=n2)
  return(out)
}
