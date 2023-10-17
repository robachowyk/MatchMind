logPossibleConfig = function(Brecords,sumD)
{
  return = 0
  if(sumD>0)
    return = sum( log(Brecords:(Brecords-sumD+1))  )
  return 
}

loglik = function(LLL, LLA, LLB, D, links, sumRowD, sumColD, gamma)
{
  # Sanity check for using logPossibleConfig(.) below
  if(ncol(D) - sum(D) + 1 <= 0){
    print("The number of records in B has to be >= to the number of linked records.")
    print("Number of records in B:")
    print(ncol(D))
    print("Number of linked records:")
    print(sum(D))
    print("The problem may come from the fact that file A is bigger than file B, which should not happen.")
    print("Number of records in A:")
    print(nrow(D))
  }
  logPossD = sum(log(gamma) * sumRowD + log(1-gamma) * (1-sumRowD)) - logPossibleConfig(ncol(D),sum(D))
  logPossD + sum(LLA[sumRowD==0]) + sum(LLB[sumColD==0]) + sum(LLL[links])
}

simulateH = function(data, D, links, omegaData, sumRowD, sumColD, eta, omega, phi)
{
  nonlinkedA = sumRowD==0
  nonlinkedB = sumColD==0
  truePIVs = sampleH(nA=dim(data$pivsA), nB=dim(data$pivsB), links=links, omegaData=omegaData, pivs_stable=data$pivs_stable, pivsA=data$pivsA, pivsB=data$pivsB, nvalues=data$nvalues, D=D, nonlinkedA=nonlinkedA, nonlinkedB=nonlinkedB, eta=eta, omega=omega, phi=phi)
  list(truepivsA=truePIVs$truepivsA, truepivsB=truePIVs$truepivsB)
}

simulateD = function(data, D, links, sumRowD, sumColD, truepivsA, truepivsB, gamma, eta, omega, phi)
{
  # Determine which observation pairs can be linked (or not) based on the true values of the PIVS
  # -> focus on stable PIVs only this should only (that we can select with the mask)
  HA = do.call(paste, c(as.data.frame(truepivsA[,data$pivs_stable]), list(sep="_")))
  HB = do.call(paste, c(as.data.frame(truepivsB[,data$pivs_stable]), list(sep="_")))
  valuesH = unique(c(HA,HB))
  HA = as.numeric(factor(HA,levels=valuesH))
  HB = as.numeric(factor(HB,levels=valuesH))
  select = F1(HA, HB, length(valuesH))
  select = F11(select, length(select))
  select = do.call("rbind", select)
  pLink = rep(gamma, nrow(data$pivsA))
  # What should we add/substract from the loglikelihood if an observation is a not linked
  LLA = rep(0,nrow(data$pivsA))
  for(k in 1:length(data$nvalues))
  {
    logpTrue = log(eta[[k]])[truepivsA[,k]]
    pMissingA = phi[[k]][2]
    pTypoA = (1-pMissingA) * (1-phi[[k]][1]) / (data$nvalues[k]-1)
    pAgreeA = (1-pMissingA) * phi[[k]][1]
    # Contribution to the likelihood
    contr = rep(pAgreeA, nrow(data$pivsA))
    contr[data$pivsA[,k] != truepivsA[,k]] = pTypoA 
    contr[data$pivsA[,k] == 0] = pMissingA
    LLA = LLA + logpTrue + log(contr)
  } 
  # What should we add/substract from the loglikelihood if an observation is a not linked
  LLB = rep(0,nrow(data$pivsB)) 
  for(k in 1:length(data$nvalues))
  {
    logpTrue = log(eta[[k]])[truepivsB[,k]]
    pMissingB = phi[[k]][3]
    pTypoB = (1-pMissingB) * (1-phi[[k]][1]) / (data$nvalues[k]-1)
    pAgreeB = (1-pMissingB) * phi[[k]][1]
    # Contribution to the likelihood
    contr = rep(pAgreeB, nrow(data$pivsB))
    contr[data$pivsB[,k] != truepivsB[,k]] = pTypoB
    contr[data$pivsB[,k] == 0] = pMissingB
    LLB = LLB + logpTrue + log(contr)
  } 
  # What do we add/substract from the loglikelihood if a pair is linked
  LLL = Matrix(0, nrow=nrow(data$pivsA), ncol=nrow(data$pivsB))
  for(k in 1:length(data$nvalues))
  {
    HA = truepivsA[ select[,1],k ]
    HB = truepivsB[ select[,2],k ]
    logpTrue = log(eta[[k]])[HA]
    pMissingA = phi[[k]][2]
    pTypoA = (1-pMissingA) * (1-phi[[k]][1]) / (data$nvalues[k]-1)
    pAgreeA = (1-pMissingA) * phi[[k]][1]
    pMissingB = phi[[k]][3]
    pTypoB = (1-pMissingB) * (1-phi[[k]][1]) / (data$nvalues[k]-1)
    pAgreeB = (1-pMissingB) * phi[[k]][1]
    # Contribution to the likelihood of linked observation from A
    helpA = rep(pAgreeA, length(HA))
    helpA[data$pivsA[select[,1],k] != HA] = pTypoA
    helpA[data$pivsA[select[,1],k] == 0] = pMissingA
    # Contribution to the likelihood of linked observation from B
    helpB = rep(pAgreeB, length(HB))
    helpB[data$pivsB[select[,2],k] != HB] = pTypoB
    helpB[data$pivsB[select[,2],k] == 0] = pMissingB
    LLL[select] = LLL[select] + logpTrue + log(helpA) + log(helpB)
    # Add unstable part if unstable
    if(!data$pivs_stable[k])
    {
      times = abs(data$regTimeB[select[,2]] - data$regTimeA[select[,1]])
      Xomega = cbind( times, 
                      data$XA[select[,1], data$pSameH.varA[[k]], drop=FALSE], 
                      data$XB[select[,2], data$pSameH.varB[[k]], drop=FALSE] )
      pSameH = exp(- as.matrix(Xomega) %*% omega[[k]])
      helpH = pSameH^(HA==HB) * ((1-pSameH)/(data$nvalues[k]-1))^(HA!=HB)
      LLL[select] = LLL[select] + log(helpH)
    }
  }  
  LLL[select][is.na(LLL[select])] = Inf
  # Complete data likelihood
  LL0 = loglik(LLL=LLL, LLA=LLA, LLB=LLB, D=D, links=links, sumRowD=sumRowD, sumColD=sumColD, gamma=gamma)
  # Single run through D
  Dsample = sampleD(S=as.matrix(select),
                    LLA=LLA, 
                    LLB=LLB, 
                    LLL=LLL,
                    gamma=pLink, 
                    loglik=LL0, 
                    D=D, 
                    nlinkrec=as.integer(sum(D)), 
                    sumRowD=sumRowD>0, 
                    sumColD=sumColD>0)
  # Sanity check: does it give the same likelihood?
  if (round(Dsample$loglik, digits = 3) != round(loglik( LLL=LLL, LLA=LLA, LLB=LLB, D=Dsample$D, links=Dsample$links+1, sumRowD=Dsample$sumRowD, sumColD=Dsample$sumColD, gamma=pLink), digits = 3))
  {
    print( "Sanity check failed" )
    print( "Log likelihood associated with new Delta:" )
    print( Dsample$loglik )
    print( "Log likelihood computed on the new Delta:" )
    print( loglik( LLL=LLL, LLA=LLA, LLB=LLB, D=Dsample$D, links=Dsample$links+1, sumRowD=Dsample$sumRowD, sumColD=Dsample$sumColD, gamma=pLink) )
  }
  Dsample
}

stEM = function(data, nIter, nBurnin, MStepIter, trace=1)
{
  for(k in 1:length(data$nvalues)){
    if(length(data$pSameH.varA[[k]])>0)
      assert( "Some variables from A to include in the survival model for unstable PIVs do not exist.", pSameH.varA[[k]] %in% colnames(XA) )
    if(length(data$pSameH.varB[[k]])>0)
      assert( "Some variables from B to include in the survival model for unstable PIVs do not exist.", pSameH.varB[[k]] %in% colnames(XB) )
  }
  
  # Parameters for PIVs
  
  # Parameter for the probability for a pair to be linked
  gamma = 0.99
  
  # Parameters for the distributions of the true values
  eta = lapply(data$nvalues, function(x) rep(1/x,x))
  
  # Parameters for the survival model describing pSameH for unstable PIVs (over time and potentially more covariates)
  # Number of coefficients: covariates from A + covariates from B + time
  ncoef = lapply(seq_along(data$pivs_stable), function(idx) if(data$pivs_stable[idx]){ 0 }else{ ncol(data$XA[, data$pSameH.varA[[idx]], drop=FALSE]) + ncol(data$XB[, data$pSameH.varB[[idx]], drop=FALSE]) + 1 } )
  omega = lapply(seq_along(data$pivs_stable), function(idx) if(data$pivs_stable[idx]){ c(-Inf) }else{ rep(0.01, ncoef[idx]) })
  
  # Parameters for the registration errors (agreement, missing in A, missing in B)
  phi = lapply(data$nvalues, function(x)  c(0.95,0.01,0.01))
  
  NmissingA = lapply(seq_along(data$nvalues), function(k) sum(data$pivsA[,k]==0))
  NmissingB = lapply(seq_along(data$nvalues), function(k) sum(data$pivsB[,k]==0))
  
  gamma.iter = array(NA, c(nIter, length(gamma)))
  eta.iter = lapply(data$nvalues, function(x) array(NA, c(nIter, x)))
  phi.iter = lapply(data$nvalues, function(x) array(NA, c(nIter, 3)))
  omega.iter = lapply(ncoef, function(x) array(NA, c(nIter, x)))
  
  time.iter=c()
  if(trace==1)
    pb = progress_bar$new(format = "Running StEM algorithm [:bar] :percent in :elapsed",       total = nIter, clear = FALSE, width= 60)
  
  if(trace==2)
    cat("-----------------------\n")
  
  iter = 1
  
  # MC-EM iteration
  for(iter in 1:nIter)
  {
    tijdM = Sys.time()
    
    D = Matrix(0, nrow=nrow(data$pivsA), ncol=nrow(data$pivsB), sparse=TRUE)
    links = which(D==1, arr.ind=TRUE)
    sumRowD = rowSums(D)
    sumColD = colSums(D)
    omegaData = list()
    
    # Burn-in period Gibbs sampler
    for(j in 1:nBurnin)
    {
      newTruePivs = simulateH(data=data, D=D, links=links, omegaData=omegaData, sumRowD=sumRowD, sumColD=sumColD, eta=eta, omega=omega, phi=phi)
      truepivsA = newTruePivs$truepivsA
      truepivsB = newTruePivs$truepivsB
      Dsample = simulateD(data=data, D=D, links=links+1, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, omega=omega, phi=phi)
      
      D = Dsample$D
      links = Dsample$links
      sumRowD = Dsample$sumRowD
      sumColD = Dsample$sumColD
      loglikelihood = Dsample$loglik
      nlinkrec = Dsample$nlinkrec
      times = abs(data$regTimeB[links[,2]+1] - data$regTimeA[links[,1]+1])
      omegaData = list()
      for(k in 1:length(data$nvalues)){
        if(nrow(links)>0){
          if(!data$pivs_stable[k]){
            Xomega = cbind( times, 
                            data$XA[links[,1]+1, data$pSameH.varA[[k]], drop=FALSE], 
                            data$XB[links[,2]+1, data$pSameH.varB[[k]], drop=FALSE] )
            omegaData[[k]] = as.matrix(Xomega) %*% omega[[k]]}}}
    }
    
    # Administration for the M-step
    Vgamma = c()
    Veta = lapply(data$nvalues, function(x) c())
    createOmegadata <- function(ncoef, stable){
      if (!stable){ return ( data.frame(matrix(nrow = 0, ncol = ncoef)) ) } }
    Vomega = mapply(createOmegadata, ncoef, data$pivs_stable, SIMPLIFY=FALSE)
    Vphi = lapply(data$nvalues, function(x) c())
    
    for(j in 1:MStepIter)
    {
      newTruePivs = simulateHchronoCPP(data=data, D=D, links=links, omegaData=omegaData, sumRowD=sumRowD, sumColD=sumColD, eta=eta, omega=omega, phi=phi)
      truepivsA = newTruePivs$truepivsA
      truepivsB = newTruePivs$truepivsB
      Dsample = simulateD(data=data, D=D, links=links+1, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, omega=omega, phi=phi)
      
      D = Dsample$D
      links = Dsample$links
      sumRowD = Dsample$sumRowD
      sumColD = Dsample$sumColD
      loglikelihood = Dsample$loglik
      nlinkrec = Dsample$nlinkrec
      times = abs(data$regTimeB[links[,2]+1] - data$regTimeA[links[,1]+1])
      omegaData = list()
      for(k in 1:length(data$nvalues)){
        if(nrow(links)>0){
          if(!data$pivs_stable[k]){
            equal = truepivsA[links[,1]+1,k] == truepivsB[links[,2]+1,k]
            Xomega = cbind( times, 
                            data$XA[links[,1]+1, data$pSameH.varA[[k]], drop=FALSE], 
                            data$XB[links[,2]+1, data$pSameH.varB[[k]], drop=FALSE] )
            # Update omega
            Vtmp_k = cbind( equal,
                            Xomega )
            Vomega[[k]] = rbind(Vomega[[k]], Vtmp_k)
            omegaData[[k]] = as.matrix(Xomega) %*% omega[[k]]}}}
      
      # Update gamma
      Vgamma[j] = sum(D)
      
      # Update eta
      for(k in 1:length(data$nvalues))
      {
        facpivsA = factor(truepivsA[,k], levels=1:data$nvalues[k])
        facpivsB = factor(truepivsB[,k], levels=1:data$nvalues[k])
        Veta[[k]] = rbind(Veta[[k]],
                          table(facpivsA[sumRowD==0]) + table(facpivsB[sumColD==0]) + table(facpivsA[sumRowD==1]))
      }

      # Update phi: count agreements / missings
      for(k in 1:length(data$nvalues))
        Vphi[[k]] = rbind( Vphi[[k]], 
                           c( sum(truepivsA[,k]==data$pivsA[,k]) + sum(truepivsB[,k]==data$pivsB[,k]), 
                              NmissingA[[k]],
                              NmissingB[[k]]) )
    }
    
    # Calculate new parameters gamma/eta/phi/omega 
    
    # New gamma
    gamma = sum(Vgamma) / (MStepIter*nrow(truepivsA))
    
    # New eta
    for(k in 1:length(data$nvalues))
      eta[[k]] = colSums(Veta[[k]])/sum(Veta[[k]])
    
    # New omega
    if(nrow(links)>0){
      loglikOmega = function(logOmega, dataOmega)
      {
        S = exp( - as.matrix(dataOmega[,!names(dataOmega) %in% c("equal")]) %*% exp(logOmega) )
        -sum(log(S) * dataOmega$equal + log(1-S) * !dataOmega$equal)
      }  
      for(k in 1:length(data$nvalues))
      {
        if(!data$pivs_stable[k])
        {
          omegaInit = rep(0.01, ncoef[k])
          omega[[k]] = exp(nlminb(omegaInit, loglikOmega, control=list(trace=FALSE), dataOmega=Vomega[[k]])$par)
        }
      } 
    }
    
    # New phi
    for(k in 1:length(data$nvalues))
    {
      Ntotal = MStepIter * (nrow(data$pivsA) + nrow(data$pivsB))
      phi[[k]][1] = sum(Vphi[[k]][,1]) / (Ntotal - sum(Vphi[[k]][,2]) - sum(Vphi[[k]][,3]))
      phi[[k]][2] = sum(Vphi[[k]][,2]) / (MStepIter*nrow(data$pivsA))
      phi[[k]][3] = sum(Vphi[[k]][,3]) / (MStepIter*nrow(data$pivsA))
    } 
    
    # Administration
    
    gamma.iter[iter,] = gamma
    
    for(k in 1:length(data$nvalues))
    {
      eta.iter[[k]][iter,] = eta[[k]]
      omega.iter[[k]][iter,] = omega[[k]]
      phi.iter[[k]][iter,] = phi[[k]]
    }
    
    if(trace==1)
      pb$tick()
    
    if(trace==2)
    {
      cat("-------------------------------------------------------------------------\n")
      cat("Iter:", iter, "\n")
      
      time.iter = time.iter + as.numeric(difftime(Sys.time(), tijdM, units="min"))
      cat("Approx time remaining:\t", (nIter-iter) * time.iter/iter , "mins\n")
      
      cat("gamma:\t", gamma, "\n")
      cat("phi:")
      for(i in 1:length(data$nvalues))
        cat("\t\t", phi[[i]], "\n")
      
      cat("eta:")
      for(i in 1:length(data$nvalues))
        cat("\t\t",  eta[[i]], "\n") 
      
      cat("omega:")
      for(i in 1:length(data$nvalues))
        cat("\t\t",  omega[[i]], "\n") 
    }
  }
  
  if(trace==1)
    pb$terminate()
  
  Delta = matrix(0, nrow=nrow(data$pivsA), ncol=nrow(data$pivsB))
  
  gamma_avg = apply(gamma.iter, 2, function(x) mean(x))
  eta_avg = lapply(eta.iter, function(x) apply(x, 2, mean))
  omega_avg = lapply(omega.iter, function(x) apply(x, 2, mean))
  phi_avg = lapply(phi.iter, function(x) apply(x, 2, mean))
  
  for(m in 1:250)
  {
    newTruePivs = simulateH(data=data, D=D, links=links, omegaData=omegaData, sumRowD=sumRowD, sumColD=sumColD, eta=eta, omega=omega, phi=phi)
    truepivsA = newTruePivs$truepivsA
    truepivsB = newTruePivs$truepivsB
    Dsample = simulateD(data=data, D=D, links=links+1, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, omega=omega, phi=phi)
    
    D = Dsample$D
    links = Dsample$links
    sumRowD = Dsample$sumRowD
    sumColD = Dsample$sumColD
    loglikelihood = Dsample$loglik
    nlinkrec = Dsample$nlinkrec
    times = abs(data$regTimeB[links[,2]+1] - data$regTimeA[links[,1]+1])
    omegaData = list()
    for(k in 1:length(data$nvalues)){
      if(nrow(links)>0){
        if(!data$pivs_stable[k]){
          Xomega = cbind( times, 
                          data$XA[links[,1]+1, data$pSameH.varA[[k]], drop=FALSE], 
                          data$XB[links[,2]+1, data$pSameH.varB[[k]], drop=FALSE] )
          omegaData[[k]] = as.matrix(Xomega) %*% omega[[k]]}}}
    
    Delta = Delta + D
  } 
  Delta = Delta / 250
  
  browseURL('https://www.youtube.com/watch?v=NTa6Xbzfq1U')
  
  list(Delta=Delta, gamma=gamma.iter, eta=eta.iter, omega=omega.iter, phi=phi.iter)
}