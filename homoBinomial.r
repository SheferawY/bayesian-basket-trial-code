homoBinomial<-function(triple,cutFisher=0.80) 
{
  # merge groups with homogeneous binomial proportions 
  # input: 
  #   data: 2 by k matrix, first row # of successes, second row # of failures
  #   groupID: group indices of each column, staring with '1','2',...,'k'
  #   toMerge: TRUE or FALSE
  # example:
  #   data=matrix(c(11,15,17,34,36,18),2);groupID=as.list(1:3);toMerge=TRUE         
  #   triple=list(data=matrix(c(11,15,17,34,36,18),2),groupID=as.list(1:3),toMerge=TRUE)
  
  while(triple$toMerge)
  {
    prop=triple$data[1,]/(triple$data[1,]+triple$data[2,])     # proportion
    sortProp=sort(prop,index.return=TRUE) # sort 
    minIX=which.min(diff(sortProp$x))     # nearest neighbor ID
    dataFisher=triple$data[,sortProp$ix[minIX:(minIX+1)]] # extract 2*2 table                 
    pval=fisher.test(dataFisher)$p.value                  # Fisher test for 2*2 table
    if(pval>=cutFisher)   
    {
      triple$data[,sortProp$ix[minIX]]=triple$data[,sortProp$ix[minIX]]+
      triple$data[,sortProp$ix[minIX+1]]                           # combine 2*2 in one column
      triple$data=matrix(triple$data[,-sortProp$ix[minIX+1]],2)    # remove the other column
      triple$groupID[[sortProp$ix[minIX]]]=union(triple$groupID[[sortProp$ix[minIX]]],
                                                 triple$groupID[[sortProp$ix[minIX+1]]]) 
      # union group ID in one group
      triple$groupID=triple$groupID[-sortProp$ix[minIX+1]] # remove the other group
      
      if(length(triple$groupID)>1){             # still more than one group
        triple=homoBinomial(triple)}            # repeat merging
      else {triple$toMerge=FALSE}               # one group already
    }
    else{triple$toMerge=FALSE                   # pval<=0.05
    break
    }
  } # end of while
  
  return(triple) # return more information
}