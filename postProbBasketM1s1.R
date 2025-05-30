postProbBasketM1s1<-function(PI,pM,n,cutEff=0.761, aP=1, bP=1) 
{ # return PS_j, PCS and PCS-k of the proposed method for basket trial design 
  # input: PI           = proportions of cancer types in sub-pop
  #        pM           = true response rate
  #        n            = sample size for single stage (i.e. n=n1+n2 )
  #        alphaFut     = early stop error rate
  #        cutEff       = cutoff for efficacy to achieve PE=alphaEff=0.30 under pjk=q1 for all (j,k)th subgroup
  #        aP,bP        = priors for for response rate p_{jk}
  # output:PCS          = PCS of optimal cutoff
  #        MSE          = Mean square error
  #        PCSk         = PCS of optimal cutoff for each cancer type
  #        PE           = emperical individual prob of efficacy
  
  JK=dim(PI); J=JK[1]; K=JK[2]  # # of cutoffs &  # of cancer types
  q0=0.20; q1=0.35              # resp rate cutoff for futility and efficacy
  numSimu=10000                 # num. of simulation 
  
  # initialize
  postProb1=postProb2=nonFutFLAG=effFLAG=Q=QTRUE=matrix(0,J,K)                   
  rowPostProb=rowEffFLAG=rep(0,J); jStar_k=jStarHat_k=colPostProb=colEffFLAG=rep(0,K)
  jStarHat=rep(0,numSimu); jStarHat_k=matrix(0,numSimu,K)
  
  cFutMatrix=array(FALSE,dim=c(numSimu,J,K))          # collection of futility matrix
  cEffMatrix=array(FALSE,dim=c(numSimu,J,K))          # collection of efficacy matrix
  cN=array(FALSE,dim=c(numSimu,J,K))                  # collection of samples
  
  # compute true jStar and true jStar-k
  for(j in 1:J) for (k in 1:K) QTRUE[j,k]=PI[j,k]*(pM[j,k]>q1)
  if(max(QTRUE)!=0){jStar=which.max(rowSums(QTRUE))}else{jStar=NA}
  for(k in 1:K) if(max(QTRUE[,k])!=0){jStar_k[k]=which.max(QTRUE[,k])}else{jStar_k[k]=NA}
  
  for (iSimu in 1:numSimu) {            # main loop
    
    N=n*matrix(rep(1,J*K),J,K) 
    #set.seed(100+iSimu)                # the value to get the same result for checking the difference
    R=apply(pM, 1:2,rbinom,n=1,size=n)  # response matrix in singel stage 
    
    # post prob in single stage 
    for (j in 1:J) for (k in 1:K) 
    {postProb2[j,k]=pbeta(q1,aP+R[j,k],bP+N[j,k]-R[j,k],lower.tail = FALSE)}  
    effFLAG=(postProb2>cutEff)*1         # flag matrix of efficacy after stage 2
    
    # final product of pi-jk with indicator functions
    for (j in 1:J) for (k in 1:K) Q[j,k]=PI[j,k]*effFLAG[j,k]    
    
    # compute estimated jStar and jStar_k
    #if(max(Q)!=0){jStarHat[iSimu]=which.max(rowSums(Q))}else{jStarHat[iSimu]=NA}
    #1#if(max(Q)!=0){if(length(unique(rowSums(Q)))!=1){jStarHat[iSimu]=which.max(rowSums(Q))}else{jStarHat[iSimu]=3}}else{jStarHat[iSimu]=NA}  
    if(max(Q)!=0){
      if( length(unique(rowSums(Q)))!=1){
        if( length(unique(rowSums(Q)))==3){jStarHat[iSimu]=which.max(rowSums(Q))}else{jStarHat[iSimu]=3} 
        if( (length(which(rowSums(Q)==max(rowSums(Q))))==2)& (which(rowSums(Q)==max(rowSums(Q)))[2]==3) ){jStarHat[iSimu]=3}
        else{if( (length(which(rowSums(Q)==max(rowSums(Q))))==2)&(which(rowSums(Q)==max(rowSums(Q)))[2]==2) ){jStarHat[iSimu]=2}else{jStarHat[iSimu]=which.max(rowSums(Q))}}
      }else{jStarHat[iSimu]=3}
    }else{jStarHat[iSimu]=NA} 
    
    for(k in 1:K ) if(max(Q[,k])!=0){jStarHat_k[iSimu,k]=which.max(Q[,k])}else{jStarHat_k[iSimu,k]=NA}
    
    cEffMatrix[iSimu,1:J,1:K]=effFLAG
    cN[iSimu,1:J,1:K]=N
  } # end of main loop
  jstarHatExist = (numSimu-sum(is.na(jStarHat)))
  MSE=sum((jStarHat-jStar)^2, na.rm = TRUE)/jstarHatExist
  PCS=sum(jStarHat==jStar,na.rm=TRUE)/jstarHatExist # one scalar
  PCSk=colSums(jStarHat_k==matrix(jStar_k,numSimu,K,byrow = TRUE),na.rm=TRUE)/numSimu   # 1 by K vector
  PE=apply(cEffMatrix,c(2,3),mean)
  sampleSize=sum(apply(cN,c(2,3),mean))
  
  return(result=list(PCS=PCS,MSE=MSE,PCSk=PCSk,PE=PE))
}