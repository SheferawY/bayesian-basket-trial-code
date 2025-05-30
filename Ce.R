Ce=function(PI,n1=20, n2=20, q0=0.20, q1=0.35,alphaFut=0.20, alphaEff=0.30, cutFisher=0.80, numSimu=1000000){
  # return: upper alphaFut and alphaEff quantile of the upper imcomplete beta distribution for method I and II
  # PI    : MIXED prevalence of each cancer under each cutoff
  # example: cutEff=Ce(n1=20, n2=20, q0=0.20, q1=0.35,alphaFut=0.20, alphaEff=0.30, cutFisher=0.80, numSimu=1000000)
  
  # return: upper alphaFut quantile of the upper imcomplete beta distribution 
  r=qbinom(1-alphaFut, n1, q0, lower.tail = TRUE)
  cutFut = pbeta((q0+q1)/2, 1+r, 1+n1-r, lower.tail = FALSE)
  
  if(alphaFut==1){cutFut=1}                                     # Single-stage n=n1    by prof Xu cutFut formula
  if(alphaFut==0){cutFut=-0.001}  #-0.001 is to enclude even ze # Single-stage n=n1+n2 by prof Xu cutFut formula
  
  JK=dim(PI); J=JK[1]; K=JK[2]  # # of cutoffs &  # of cancer types
  q0=0.20; q1=0.35                    # resp rate cutoff for futility and efficacy
  
  #pM0=matrix(rep(q0,J*K),J,K)        # true p_jk=q1 for all (j,k) cells under scenario 1
  pM=matrix(rep(q1,J*K),J,K)          # true p_jk=q1 for all (j,k) cells under scenario 1
  
  # approximation of beta priors for marginal responses pj. and p.k
  aj.=rep(0,J); b.k=rep(0,K)
  for (j in 1:J){ for (k in 1:K){ 
    aj.[j]={(3*sum(PI[j,])^2)/(2*sum(PI[j,]^2))}-0.5 
    b.k[k]={(3*sum(PI[,k])^2)/(2*sum(PI[,k]^2))}-0.5
  }}
  
  # initialize:create array and collection of matrix
  CeI=CeII=matrix(FALSE,J,K); rowCeII=rep(FALSE,J); colCeII=rep(FALSE,K)
  upIBetaIIjk=upIBetaII=R=N=array(FALSE,dim=c(numSimu,J,K))                        
  rowR=rowN=rowPostProb=rowEffFLAG=matrix(FALSE,numSimu,J) 
  colR=colN=colPostProb=colEffFLAG=matrix(FALSE,numSimu,K)
  
  R1=apply(pM, 1:2,rbinom,n=numSimu,size=n1)  # response matrix in stage 1
  R2=apply(pM, 1:2,rbinom,n=numSimu,size=n2)  # response matrix in stage 2
  delta = pbeta((q0+q1)/2, 1+R1, 1+n1-R1, lower.tail = FALSE)>cutFut    # an array of promising indicator
  
  for (i in 1:numSimu){
    # post prob after stage 2
    for (j in 1:J){ for (k in 1:K){ 
      R[i,j,k]=R1[i,j,k]+R2[i,j,k]*delta[i,j,k]   # total responses matrix after stage 2
      N[i,j,k]=n1+n2*delta[i,j,k]                 # sample size matrix after stage 2
      upIBetaII[i,j,k] = pbeta(q1, 1+R[i,j,k], 1+N[i,j,k]-R[i,j,k], lower.tail = FALSE)
    }}
    
    # marginal responses and sample sizes 
    colR[i,]=colSums(R[i,,]); rowR[i,]=rowSums(R[i,,]); colN[i,]=colSums(N[i,,]); rowN[i,]=rowSums(N[i,,])
    
    # marginal post prob by row   [most tricky code here!!]
    triple=list(data=rbind(rowR[i,],rowN[i,]-rowR[i,]),groupID=as.list(1:J),toMerge=TRUE)   
    res=list(data=rbind(rowR[i,],rowN[i,]-rowR[i,]),groupID=as.list(1:J),toMerge=TRUE)#=homoBinomial(triple,cutFisher)
    for (iG in 1:length(res$groupID))     # row margin post prob # Approximating beta prior distribution for pj. and p.k
    { rowPostProb[i,res$groupID[[iG]]]=pbeta(q1,aj.[res$groupID[[iG]]]+res$data[1,iG],    
                                             aj.[res$groupID[[iG]]]+res$data[2,iG],
                                             lower.tail = FALSE)    }
    
    # marginal post prob by column
    triple=list(data=rbind(colR[i,],colN[i,]-colR[i,]),groupID=as.list(1:K),toMerge=TRUE)
    res=list(data=rbind(colR[i,],colN[i,]-colR[i,]),groupID=as.list(1:K),toMerge=TRUE)#=homoBinomial(triple,cutFisher)
    for (iG in 1:length(res$groupID))     # row margin post prob  # Approximating beta prior distribution for pj. and p.k
    { colPostProb[i,res$groupID[[iG]]]=pbeta(q1,b.k[res$groupID[[iG]]]+res$data[1,iG],  
                                             b.k[res$groupID[[iG]]]+res$data[2,iG],
                                             lower.tail = FALSE)    }
    
    for (j in 1:J){
      for (k in 1:K){
        upIBetaIIjk[i,j,k]=upIBetaII[i,j,k]*apply(rbind(rowPostProb[i,j],colPostProb[i,k]),2,max)
      }}
  }
  
  # find the matrix of threshold for method I and II
  for (j in 1:J){
    for (k in 1:K){
      rowCeII[j] = quantile(rowPostProb[,j], probs=1-alphaEff, names = FALSE)       # row Ce for method II two-stage
      colCeII[k] = quantile(colPostProb[,k], probs=1-alphaEff, names = FALSE)       # column Ce for method II two-stage
      CeI[j,k]   = quantile(upIBetaII[,j,k], probs=1-alphaEff, names = FALSE)       # Ce for method I two-stage
      CeII[j,k]  = quantile(upIBetaIIjk[,j,k], probs=1-alphaEff, names = FALSE)     # Ce for method II two-stage
    }} 
  
  return(list(cutFut=cutFut,
              CeI=CeI,
              rowCeII=rowCeII,
              colCeII=colCeII,
              CeII=CeII))
}