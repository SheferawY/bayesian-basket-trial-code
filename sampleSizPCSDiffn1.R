sampleSizPCSDiffn1<-function(PI,pM,alphaFut=0.20,cutEff_In1,cutEff_In1II,row_cutEff1,col_cutEff1) 
{ # cutEff_In1 is CeI for method I, cutEff_In1II is CeI for method II (it's equal to cutEff_In1, when we approx to make equal alphaEff)
  # ,row_cutEff1 and col_cutEff1 are vectors 
  # return the PCS and total sample size for method I and II at different n1 and fixd n2=20
  # example: s1_pcsDiffn1=sampleSizPCSDiffn1(PI,pM,alphaFut=0.20,cutEff_In1,row_cutEff1,col_cutEff1)
  n1=seq(10,30,1); n2=20   
  PCSM1=PCSM2=cutEff_I=meanPE_I=meanPE_II=rep(0,length(n1)) 
  for (i in 1:length(n1)) { # length(n1)
    out_M1=postProbBasketM1(PI,pM,n1=n1[i],n2,alphaFut,cutEff=cutEff_In1[i])        #  method I at any scen
    out_M2=postProbBasketM2(PI,pM,n1=n1[i],n2,alphaFut,cutEff=cutEff_In1II[i],
                            row_cutEff=row_cutEff1[i,],col_cutEff=col_cutEff1[i,])          #  method II at any scen 
    
    cutEff_I[i]=out_M1$cutEff  # cutEff for method I at any scen
    PCSM1[i]=out_M1$PCSd       # PCS for method I at any scen
    PCSM2[i]=out_M2$PCSd       # PCS for method II at any scen
    meanPE_I[i]=mean(out_M1$PE)  # mean PE cutEff for method I used to check equality of mean PE for Method and I under scen 1
    meanPE_II[i]=mean(out_M2$PE)  # mean PE cutEff for method I used to check equality of mean PE for Method and I under scen 1
  }
  out<-list(n1=n1,n2=n2,cutEff_I=cutEff_I,PCSM1=PCSM1,PCSM2=PCSM2,meanPE_I=meanPE_I,meanPE_II=meanPE_II)
  return(out)
}
