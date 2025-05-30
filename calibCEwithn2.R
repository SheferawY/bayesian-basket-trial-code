calibCEwithn2<-function(PI,n1=20,q0=0.20, q1=0.35,alphaFut=0.20,alphaEff=0.30,cutFisher=0.80,numSimu=1000000) 
{ # generate Ce at fixed n2=c(10:30) and fixed n1=20 UNDER TWO-STAGE
  n2=seq(10,30,1)
  JK=dim(PI); J=JK[1]; K=JK[2]
  Cf_n2=rep(0,length(n2)) 
  cutEff_In2=array(0,dim=c(length(n2),J,K));  cutEff_IIn2=array(0,dim=c(length(n2),J,K))                        
  row_cutEff=matrix(0,length(n2),J);   col_cutEff=matrix(0,length(n2),K)
  for (i in 1:length(n2)) { # length(n2)
    
     Ce_n2i=CeIIprof(n1,n2=n2[i],q0, q1,alphaFut,alphaEff,cutFisher,numSimu) #  all methods calibrated Cf and Ce for each n2 and fixed n1
    # threshold outputs
    Cf_n2[i]=Ce_n2i$cutFut          # calibrated Cf for each n2 at fixed n1
    cutEff_In2[i,,]=Ce_n2i$CeI      # calibrated CeI for each n2 at fixed n1
    row_cutEff[i,]=Ce_n2i$rowCeII   # calibrated rowCeII for each n2 at fixed n1
    col_cutEff[i,]=Ce_n2i$colCeII   # calibrated colCeII for each n2 at fixed n1
    cutEff_IIn2[i,,]=Ce_n2i$CeII    # calibrated CeII for each n2 at fixed n1
  }
  out<-list(Cf_n2=Cf_n2,cutEff_In2=cutEff_In2,
            row_cutEff=row_cutEff,col_cutEff=col_cutEff,
            cutEff_IIn2=cutEff_IIn2,
            n1=n1,n2=n2)
  return(out)
}