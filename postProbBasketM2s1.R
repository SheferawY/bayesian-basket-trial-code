postProbBasketM2s1<-function(PI,pM,n,cutEff=0.761,
                             row_cutEff=c(0.7984352, 0.8025440, 0.7889145),
                             col_cutEff=c(0.7796536, 0.7759109, 0.7804498, 0.7759109),
                             cutFisher=0.80, aP=1, bP=1)
{ # input: PI           = proportions of cancer types in sub-pop
  #        pM           = true response rate
  #        n            = sample size for single stage (i.e. n=n1+n2)
  #        alphaFut     = early stop error rate
  #        cutEff       = cutoff for efficacy to achieve PE=alphaEff=0.30 under pjk=q1 for all (j,k)th subgroup
  #        row_cutEff   = row marginal cutoff for efficacy to achieve PE=alphaEff=0.30 under pjk=q1 for all (j,k)th subgroup
  #        col_cutEff   = column marginal cutoff for efficacy to achieve PE=alphaEff=0.30 under pjk=q1 for all (j,k)th subgroup
  #        cutFisher    = cutoff for Fisher test
  #        aP,bP        = priors for for response rate p_{jk}
  # output:PCS          = PCS of optimal cutoff
  #        MSE          = Mean square error
  #        PCSk         = PCS of optimal cutoff for each cancer type
  #        PE           = emperical individual prob of efficacy
  
  JK=dim(PI); J=JK[1]; K=JK[2]  # # of cutoffs &  # of cancer types
  q0=0.20; q1=0.35              # resp rate cutoff for futility and efficacy
  numSimu=10000                 # num. of simulation
  
  # approximation of beta priors for marginal responses pj. and p.k
  var_pMj=aj.=rep(0, J) # variance and priors for jth row
  var_pMk=b.k=rep(0, K) # variance and priors for kth column
  for (j in 1:J){ for(k in 1:K){
    var_pMj[j]=(sum(PI[j,]^2))*(aP*bP)/((aP+bP)^2*(aP+bP+1)*sum(PI[j,])^2 )
    aj.[j]=(1-4*var_pMj[j])/(8*var_pMj[j])
    var_pMk[k]=(sum(PI[,k]^2))*(aP*bP)/((aP+bP)^2*(aP+bP+1)*sum(PI[,k])^2 )
    b.k[k]=(1-4*var_pMk[k])/(8*var_pMk[k])
  }}
  
  # initialize
  postProb1=postProb2=nonFutFLAG=effFLAG=Q=QTRUE=effFLAG_enhanced=matrix(0,J,K)    
  rowPostProb=rowEffFLAG=rep(0,J); jStar_k=jStarHat_k=colPostProb=colEffFLAG=rep(0,K)
  jStarHat=rep(0,numSimu); jStarHat_k=matrix(0,numSimu,K)
  
  cFutMatrix=array(FALSE,dim=c(numSimu,J,K))          # collection of futility matrix
  cEffMatrix=array(FALSE,dim=c(numSimu,J,K))          # collection of efficacy matrix
  cN=array(FALSE,dim=c(numSimu,J,K))                  # collection of samples
  
  ## compute true jStar and true jStar-k
  #rowPM=rowSums(pM*PI)/rowSums(PI)   # marginal response rate by row     
  #colPM=colSums(pM*PI)/colSums(PI)   # marginal response rate by column 
  
  # compute true jStar and true jStar-k
  for (j in 1:J) for (k in 1:K) QTRUE[j,k]=PI[j,k]*(pM[j,k]>q1)#*max(rowPM[j]>q1,colPM[k]>q1) 
  if(max(QTRUE)!=0){jStar=which.max(rowSums(QTRUE))}else{jStar=NA}
  for(k in 1:K) if(max(QTRUE[,k])!=0){jStar_k[k]=which.max(rowSums(PI)*pM[,k])}else{jStar_k[k]=NA}   
  
  for (iSimu in 1:numSimu) {       # main loop
    
    N=n*matrix(rep(1,J*K),J,K)
    
    R=apply(pM, 1:2,rbinom,n=1,size=n)  # response matrix in singel stage 
    
    # post prob in single stage 
    for (j in 1:J) for (k in 1:K) 
    {postProb2[j,k]=pbeta(q1,aP+R[j,k],bP+N[j,k]-R[j,k],lower.tail = FALSE)}  
    effFLAG=(postProb2>cutEff)*1          # flag matrix of efficacy after stage 2
    
    # marginal responses and sample sizes #OLD
    colR=colSums(R); rowR=rowSums(R); colN=colSums(N); rowN=rowSums(N)
    
    # marginal post prob by row   [most tricky code here!!]
    triple=list(data=rbind(rowR,rowN-rowR),groupID=as.list(1:J),toMerge=TRUE)   
    res=homoBinomial(triple,cutFisher)
    for (iG in 1:length(res$groupID))     # row margin post prob # Approximating beta prior distribution for pj.
    { rowPostProb[res$groupID[[iG]]]=pbeta(q1,aj.[res$groupID[[iG]]]+res$data[1,iG],    
                                           aj.[res$groupID[[iG]]]+res$data[2,iG],
                                           lower.tail = FALSE)    }
    
    # marginal post prob by column
    triple=list(data=rbind(colR,colN-colR),groupID=as.list(1:K),toMerge=TRUE)
    res=homoBinomial(triple,cutFisher)
    for (iG in 1:length(res$groupID))     # row margin post prob # Approximating beta prior distribution for p.k
    { colPostProb[res$groupID[[iG]]]=pbeta(q1,b.k[res$groupID[[iG]]]+res$data[1,iG],  
                                           b.k[res$groupID[[iG]]]+res$data[2,iG],
                                           lower.tail = FALSE)    }
    
    rowEffFLAG=(rowPostProb>row_cutEff)*1     # flag vector of row efficacy after stage 1
    colEffFLAG=(colPostProb>col_cutEff)*1     # flag vector of column efficacy after stage 1
    
    # final product of pi-jk with indicator functions
    for (j in 1:J){
      for (k in 1:K){ 
        effFLAG_enhanced[j,k]=effFLAG[j,k]*max(rowEffFLAG[j],colEffFLAG[k]) 
        Q[j,k]=PI[j,k]*effFLAG_enhanced[j,k]
      }
    }
    
    # compute estimated jStar and jStar_k
    #if(max(Q)!=0){jStarHat[iSimu]=which.max(rowSums(Q))}else{jStarHat[iSimu]=NA}
    #1#if(max(Q)!=0){if(length(unique(rowSums(Q)))!=1){jStarHat[iSimu]=which.max(rowSums(Q))}else{jStarHat[iSimu]=3}}else{jStarHat[iSimu]=NA} # we declare j*=3 when all are the same except all (j,k) of Q are zero 
    if(max(Q)!=0){
      if( length(unique(rowSums(Q)))!=1){
        if( length(unique(rowSums(Q)))==3){jStarHat[iSimu]=which.max(rowSums(Q))}else{jStarHat[iSimu]=3} 
        if( (length(which(rowSums(Q)==max(rowSums(Q))))==2)& (which(rowSums(Q)==max(rowSums(Q)))[2]==3) ){jStarHat[iSimu]=3}
        else{if( (length(which(rowSums(Q)==max(rowSums(Q))))==2)&(which(rowSums(Q)==max(rowSums(Q)))[2]==2) ){jStarHat[iSimu]=2}else{jStarHat[iSimu]=which.max(rowSums(Q))}}
      }else{jStarHat[iSimu]=3}
    }else{jStarHat[iSimu]=NA} # we declare j*=3 when all are the same except all (j,k) of Q are zero
    
    for(k in 1:K ) if(max(Q[,k])!=0){jStarHat_k[iSimu,k]=which.max(rowSums(PI)*effFLAG_enhanced[,k])}else{jStarHat_k[iSimu,k]=NA}  # for PCSk 
    
    cEffMatrix[iSimu,1:J,1:K]=effFLAG_enhanced
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