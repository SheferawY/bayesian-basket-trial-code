Cf = function(n1=20,q0=0.20,q1=0.35,alphaFut=0.20){
  # return: upper alphaFut quantile of the upper imcomplete beta distribution 
  # n1=20;q0=0.20;q1=0.35;alphaFut=0.20
  # example: cutFut=Cf(n1=20,q0=0.20,q1=0.35,alphaFut=0.20)
  
  r=qbinom(1-alphaFut, n1, q0, lower.tail = TRUE)
  upIBeta = pbeta((q0+q1)/2, 1+r, 1+n1-r, lower.tail = FALSE)
  return(upIBeta)
  
}


#Ce UNDER SINGLE STAGE 
#r=qbinom(1-alphaEff, n1+n2, q1, lower.tail = TRUE)
#cutEff=pbeta(q1, 1+r, 1+40-r, lower.tail = FALSE)