sampleSizPCSDiffn2<-function(PI,pM,alphaFut=0.20,cutEff_In2,cutEff_In2II,row_cutEff2,col_cutEff2) 
{ # cutEff_In2 is CeI for method I, cutEff_In2II is CeI for method II ((it's equal to cutEff_In2, when we approx to make equal alphaEff))
  # , row_cutEff2 and col_cutEff2 are vectors 
  # return the PCS and total sample size for method I and II at different n2 and fixd n1=20
  # example: s1_pcsDiffn2=sampleSizPCSDiffn2(PI,pM,alphaFut=0.20,cutEff_In2,row_cutEff2,col_cutEff2)
  n2=seq(10,30,1); n1=20   
  PCSM1=PCSM2=cutEff_I=meanPE_I=meanPE_II=rep(0,length(n2)) 
  for (i in 1:length(n2)) { # length(n2)
    out_M1=postProbBasketM1(PI,pM,n1,n2=n2[i],alphaFut,cutEff=cutEff_In2[i])        #  method I at any scen
    out_M2=postProbBasketM2(PI,pM,n1,n2=n2[i],alphaFut,cutEff=cutEff_In2II[i],
                            row_cutEff=row_cutEff2[i,],col_cutEff=col_cutEff2[i,])          #  method II at any scen 
    
    cutEff_I[i]=out_M1$cutEff  # cutEff for method I at any scen
    PCSM1[i]=out_M1$PCSd       # PCS for method I at any scen
    PCSM2[i]=out_M2$PCSd       # PCS for method II at any scen
    meanPE_I[i]=mean(out_M1$PE)  # mean PE cutEff for method I used to check equality of mean PE for Method and I under scen 1
    meanPE_II[i]=mean(out_M2$PE)  # mean PE cutEff for method I used to check equality of mean PE for Method and I under scen 1
  }
  out<-list(n1=n1,n2=n2,cutEff_I=cutEff_I,PCSM1=PCSM1,PCSM2=PCSM2,meanPE_I=meanPE_I,meanPE_II=meanPE_II)
  return(out)
}
