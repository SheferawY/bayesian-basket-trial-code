postProbBasketM1<-function(PI,pM,n1,n2,alphaFut=0.20,cutEff=0.652, aP=1, bP=1)   
{ # return PS_j, PCS and PCS-k of the proposed method for basket trial design 
  # input: PI           = proportions of cancer types in sub-pop
  #        pM           = true response rate
  #        n1           = stage 1 sample size for each subgroup
  #        n2           = stage 2 sample size for non-furtile subgroup
  #        alphaFut     = early stop error rate
  #        cutEff       = cutoff for efficacy to achieve PE=alphaEff=0.30 under pjk=q1 for all (j,k)th subgroup
  #        aP,bP        = priors for for response rate p_{jk}
  # output:PCS          = PCS of optimal cutoff
  #        MSE          = Mean square error
  #        PCSk         = PCS of optimal cutoff for each cancer type
  #        PE           = emperical individual prob of efficacy
  #        R            = empirical rate of sample size saving
  #        PnonFut      = emperical individual prob of futility at stage 1
  
  JK=dim(PI); J=JK[1]; K=JK[2]        # number of cutoffs &  number of cancer types
  q0=0.20; q1=0.35                    # resp rate cutoff for futility and efficacy
  numSimu=10000                       # num. of simulation 
  
  # initialize
  postProb1=postProb2=nonFutFLAG=effFLAG=Q=QTRUE=matrix(0,J,K)                   
  rowPostProb=rowEffFLAG=rep(0,J); jStar_k=jStarHat_k=colPostProb=colEffFLAG=rep(0,K)
  jStarHat=rep(0,numSimu); jStarHat_k=matrix(0,numSimu,K)
  
  cFutMatrix=array(FALSE,dim=c(numSimu,J,K))          # collection of futility matrix
  cEffMatrix=cnonFutMatrix=array(FALSE,dim=c(numSimu,J,K))          # collection of efficacy matrix
  cN=array(FALSE,dim=c(numSimu,J,K))                  # collection of samples
  
  # calibration of cutFut
  qFut=qbinom(1-alphaFut,n1,q0)   # use instead of alpha correction
  cutFut=pbeta((q0+q1)/2,aP+qFut,bP+n1-qFut,lower.tail = FALSE)
  cutFut=0.50
  
  # compute true jStar and true jStar-k
  for(j in 1:J) for (k in 1:K) QTRUE[j,k]=PI[j,k]*(pM[j,k]>q1)
  if(max(QTRUE)!=0){jStar=which.max(rowSums(QTRUE))}else{jStar=NA}
  for(k in 1:K) if(max(QTRUE[,k])!=0){jStar_k[k]=which.max(QTRUE[,k])}else{jStar_k[k]=NA}
  
  for (iSimu in 1:numSimu) {          # main loop
    
    #set.seed(100+iSimu)                  # the value to get the same result for checking the difference
    R1=apply(pM, 1:2,rbinom,n=1,size=n1)  # response matrix in stage 1
    for (j in 1:J) for (k in 1:K)         # post prob after stage 1
    {postProb1[j,k]=pbeta((q0+q1)/2,aP+R1[j,k],bP+n1-R1[j,k],lower.tail = FALSE)}  #PRIOR alpha=beta=12 instead of 1
    nonFutFLAG=(postProb1>cutFut)*1       # flag matrix of non-futility after stage 1
    
    R2=apply(pM, 1:2,rbinom,n=1,size=n2)  # response matrix in stage 2
    
    R=R1+R2*nonFutFLAG   # total responses matrix after stage 2
    N=n1+n2*nonFutFLAG   # sample size matrix after stage 2
    
    # post prob after stage 2
    for (j in 1:J) for (k in 1:K) 
    {postProb2[j,k]=pbeta(q1,aP+R[j,k],bP+N[j,k]-R[j,k],lower.tail = FALSE)}  
    effFLAG=(postProb2>cutEff)*1          # flag matrix of efficacy after stage 2
    
    # final product of pi-jk with indicator functions
    for (j in 1:J) for (k in 1:K) Q[j,k]=PI[j,k]*effFLAG[j,k]    
    
    # compute estimated jStar and jStar_k
    if(max(Q)!=0){
      if( length(unique(rowSums(Q)))!=1){
        if( length(unique(rowSums(Q)))==3){jStarHat[iSimu]=which.max(rowSums(Q))}else{jStarHat[iSimu]=3} 
        if( (length(which(rowSums(Q)==max(rowSums(Q))))==2)& (which(rowSums(Q)==max(rowSums(Q)))[2]==3) ){jStarHat[iSimu]=3}
        else{if( (length(which(rowSums(Q)==max(rowSums(Q))))==2)&(which(rowSums(Q)==max(rowSums(Q)))[2]==2) ){jStarHat[iSimu]=2}else{jStarHat[iSimu]=which.max(rowSums(Q))}}
      }else{jStarHat[iSimu]=3}
    }else{jStarHat[iSimu]=NA} 
    
    for(k in 1:K ) if(max(Q[,k])!=0){jStarHat_k[iSimu,k]=which.max(Q[,k])}else{jStarHat_k[iSimu,k]=NA}
    
    cnonFutMatrix[iSimu,1:J,1:K]=nonFutFLAG
    cEffMatrix[iSimu,1:J,1:K]=effFLAG
    cN[iSimu,1:J,1:K]=N
  } # end of main loop
  jstarHatExist = (numSimu-sum(is.na(jStarHat)))
  MSE=sum((jStarHat-jStar)^2, na.rm = TRUE)/jstarHatExist
  PCS=sum(jStarHat==jStar,na.rm=TRUE)/jstarHatExist  # one scalar
  PCSk=colSums(jStarHat_k==matrix(jStar_k,numSimu,K,byrow = TRUE),na.rm=TRUE)/numSimu   # 1 by K vector
  PE=apply(cEffMatrix,c(2,3),mean); PnonFut=apply(cnonFutMatrix,c(2,3),mean) # percentage of efficacy and non-futile respectively
  R=1-sum(apply(cN,c(2,3),mean))/(J*K*(n1+n2))
  
  return(result=list(PCS=PCS,MSE=MSE,PCSk=PCSk,
                     PnonFut=PnonFut,PE=PE,R=R))
}