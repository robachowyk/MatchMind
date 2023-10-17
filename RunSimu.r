################################################################################
#                                 SIMULATE DATA                                #
################################################################################

library(mvtnorm) 
library(numDeriv)  
library(progress)
library(Matrix)
library(MASS) 
library(plyr)
library(testit)
library(Metrics)

seed = round(runif(1,1,100000))

set.seed(seed) 

dataA = 250
dataB = 300
Nlinks = 100

Nval       = c(     8,      2,      5,      5,     15) 
ptypos     = c( 0.025,  0.025,  0.025,  0.025,  0.025) 
pmissing1  = c(  0.05,   0.05,   0.05,   0.05,   0.02)  
pmissing2  = c(  0.05,   0.05,   0.05,   0.05,   0.02)  

# Create two files of PIVS without overlapping units
NRecords = c(dataA, dataB)

for(i in 1:2)
{	
  dataSet=c()
  for(u in 1:length(Nval))
  {
    # Simulate different probabilities for the true values (if desired)
    xp =  exp(0 *(0:(Nval[u]-1)))
    probx = xp/sum(xp)
    
    dataSet = cbind(dataSet, sample(1:Nval[u], NRecords[i], replace=TRUE, prob=probx))
  }
  dataSet = as.data.frame(dataSet)
  names(dataSet) = c(paste("V", 1:(ncol(dataSet)), sep=""))
  assign(paste("dataSet", i, sep=""),dataSet  )
}

# Add overlapping units
match = rep(FALSE,nrow(dataSet1))
match[1:Nlinks] = TRUE
dataSet2[1:Nlinks,] = dataSet1[1:Nlinks,]

# Add missings and typos
for(x in 1:ncol(dataSet1))
{ 
  biased = as.logical(rbinom(nrow(dataSet1),1,ptypos[x]))
  if(sum(biased)>0)
    dataSet1[,x][biased] = sapply(dataSet1[,x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
}  

for(x in 1:ncol(dataSet2))
{ 
  biased = as.logical(rbinom(nrow(dataSet2),1,ptypos[x]))
  if(sum(biased)>0)
    dataSet2[,x][biased] = sapply(dataSet2[,x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
} 

for(x in 1:ncol(dataSet1))
{ 
  biased = as.logical(rbinom(nrow(dataSet1),1,pmissing1[x]))
  if(sum(biased)>0)
    dataSet1[,x][biased] = NA
}  

for(x in 1:ncol(dataSet2))
{ 
  biased = as.logical(rbinom(nrow(dataSet2),1,pmissing2[x]))
  if(sum(biased)>0)
    dataSet2[,x][biased] = NA
}  

# Add covariates (if desired)
dataSet1$gender    = runif(nrow(dataSet1), 0, 1)
dataSet1$blood     = runif(nrow(dataSet1), 0, 3)
dataSet1$siblings  = runif(nrow(dataSet1), 0, 6)
dataSet2$city      = runif(nrow(dataSet2), 0, 18)
dataSet2$parents   = runif(nrow(dataSet2), 0, 2)

XA = data.frame(dataSet1[, c("gender", "blood", "siblings")])
XB = data.frame(dataSet2[, c("city", "parents")])

# Add registration time for the survival model
regTimeA = runif(nrow(dataSet1), 0, 1)
regTimeB = runif(nrow(dataSet2), 1, 2)
time_difference = regTimeB[1:Nlinks] - regTimeA[1:Nlinks]

# Add some of the external data for the survival model
gender_data = dataSet1[1:Nlinks, "gender"]

# Create a moving process for a potential unstable variable (proba of same values decreases with time)
proba_same_H = exp( - 0.3 * (time_difference + gender_data) )

# registration time is by default taken into account in the process, add the external data wanted
pSameH.varA = list(c(), c(), c(), c(), c("gender"))
pSameH.varB = list(c(), c(), c(), c(), c())

pivs_stable = c(TRUE, TRUE, TRUE, TRUE, FALSE)

x = which(!pivs_stable)
for(i in 1:Nlinks)
{
  is_not_moving = rbinom(1, 1, proba_same_H[i])
  if(!is_not_moving)
  {
    # the value MUST be different
    dataSet2[i,x] = sample((1:Nval[x])[-c(dataSet1[i,x])], 1)
  }
}  

A = dataSet1
B = dataSet2

pivs = c("V1", "V2", "V3", "V4", "V5")

# Recode the pivs
levels_pivs = lapply(pivs, function(x) levels(factor(as.character(c(A[,x], B[,x])))))

for(i in 1:length(pivs))
{
  A[,pivs[i]] = as.numeric(factor(as.character(A[,pivs[i]]), levels=levels_pivs[[i]]))
  B[,pivs[i]] = as.numeric(factor(as.character(B[,pivs[i]]), levels=levels_pivs[[i]]))
}

# Value for missing data
A[,pivs][ is.na(A[,pivs]) ] = 0
B[,pivs][ is.na(B[,pivs]) ] = 0

nvalues = sapply(levels_pivs, length)

################################################################################
#                   LAUNCH THE ALGORITHM ON SIMULATED DATA                     #
################################################################################

getwd()
dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

dataSimu = list( pivsA        = A[, pivs], 
                 pivsB        = B[, pivs], 
                 nvalues      = nvalues,
                 pivs_stable  = pivs_stable, 
                 regTimeA     = regTimeA, 
                 regTimeB     = regTimeB,
                 XA           = XA,
                 XB           = XB,
                 pSameH.varA  = pSameH.varA,
                 pSameH.varB  = pSameH.varB )

source("recordlinkage.r")
Rcpp:::sourceCpp("functions.cpp")

fit = stEM( data       = dataSimu, 
            nIter      = 250, 
            nBurnin    = 25, 
            MStepIter  = 25, 
            trace      = 1 )

################################################################################
#                            VISUALISE THE RESULTS                             #
################################################################################

par(mfrow=c(1,2))
# True Delta underlying the data
Delta = matrix(0, nrow=nrow(dataSimu$pivsA), ncol=nrow(dataSimu$pivsB))
for (l in 1:Nlinks)
{
  Delta[l,l]=1
}
image(x=1:nrow(dataSimu$pivsB), y=1:nrow(dataSimu$pivsA), z=t(as.matrix(Delta)), xlab="obs. in B", ylab="obs. in A")
title(main = c("True linkage matrix", "underlying the data", sprintf("%s linked record pairs", sum(Delta))), font.main = 4)
# Estimated Delta 
image(x=1:nrow(dataSimu$pivsB), y=1:nrow(dataSimu$pivsA), z=t(as.matrix(fit$Delta)), xlab="obs. in B", ylab="obs. in A")
if(any(!dataSimu$pivs_stable)){
  # at least 1 PIV is considered unstable
  title(main = c("Estimated linkage matrix fitting the data","taking into account unstability in PIVs", sprintf("%s linked record pairs (proba > .5)", sum(fit$Delta>0.5))), font.main = 4)
}else{
  # all PIVs are considered stable
  title(main = c("Estimated linkage matrix fitting the data","not taking into account", "unstability in PIVs", sprintf("%s linked record pairs (proba > .5)", sum(fit$Delta>0.5))), font.main = 4)
}

par(mfrow=c(1,1))
# gamma
gamma = fit$gamma
plot(gamma, type="l", ylim=c(0,1), xlab = "MC-EM Iterations", ylab = "gamma")
abline(h = Nlinks/dataA, col="red")

par(mfrow=c(2,3))
# eta
eta = fit$eta
for(k in 1:length(eta))
{
  vec = c("eta 1")
  plot(eta[[k]][,1], ylim=c(0,1), type="l", col=1, xlab = "MC-EM Iterations", ylab = sprintf("PIV %s", pivs[k]))
  for(j in 2:ncol(eta[[k]]))
  {
    lines(eta[[k]][,j], ylim=c(0,1), col=j)
    vec <- append(vec, sprintf("eta %s", j))
  }
  legend("topright", legend=vec, col=seq_len(ncol(eta[[k]])), lty=1:4, cex=0.8)
} 

par(mfrow=c(2,3))
# omega
omega = fit$omega
for(k in 1:length(omega))
{
  if(!dataSimu$pivs_stable[k]){
    vec = c("time difference")
    plot(omega[[k]][,1], ylim=c(min(omega[[k]])-0.5, max(omega[[k]])+0.5), type="l", col=1, xlab = "MC-EM Iterations", ylab = sprintf("Survival model coef. for unstable PIV %s", pivs[k]))
    if(ncol(omega[[k]])>=2){
      vec <- append(vec, pSameH.varA[[k]])
      vec <- append(vec, pSameH.varB[[k]])
      for(c in 2:ncol(omega[[k]]))
      {
        lines(omega[[k]][,c], col=c)
      }
    }
    legend("topright", legend=vec, col=c(1:ncol(omega[[k]])), lty=1:4, cex=0.8)
  }
}

par(mfrow=c(2,3))
# phi
phi = fit$phi
for(k in 1:length(phi))
{
  plot(phi[[k]][,1], ylim=c(0,1), type="l", col=1, xlab = "MC-EM Iterations", ylab = sprintf("PIV %s", pivs[k]))
  lines(phi[[k]][,2], col=2)
  lines(phi[[k]][,3], col=3)
  legend("topright", legend=c("agree", "missings A", "missings B"), col=c(1,2,3), lty=1:4, cex=0.8)
}

# Confusion matrix
FPR_list = c()
TPR_list = c()
for(t in seq(0, 1.1, by=0.1)){
  TP_ = sum( (fit$Delta > t) & (Delta == 1) )
  TN_ = sum( (fit$Delta < t) & (Delta == 0) )
  FP_ = sum( (fit$Delta > t) & (Delta == 0) )
  FN_ = sum( (fit$Delta < t) & (Delta == 1) )
  FPR_list = append(FPR_list, FP_ / (FP_ + TN_))
  TPR_list = append(TPR_list, TP_ / (TP_ + FN_))
}
TruePositive = (fit$Delta > 0.5) & (Delta == 1)
TrueNegative = (fit$Delta < 0.5) & (Delta == 0)
FalsePositive = (fit$Delta > 0.5) & (Delta == 0)
FalseNegative = (fit$Delta < 0.5) & (Delta == 1)
TP = sum( TruePositive )
TN = sum( TrueNegative )
FP = sum( FalsePositive )
FN = sum( FalseNegative )
FPR = FP / (FP + TN)
TPR = TP / (TP + FN)
par(mfrow=c(1,2))
plot(FPR_list, TPR_list, xlim=c(0,1), ylim=c(0,1))
title(sprintf("AUC: %s", round(auc(Delta, fit$Delta),2)))
lines(c(0,1), c(0,1))
plot(FPR_list, TPR_list, xlim=c(0,0.002), ylim=c(0,1))
title(c(sprintf("TP: %s, FP: %s, FN: %s", TP, FP, FN), sprintf("FPR: %s, TPR: %s", round(FPR,2), round(TPR,2)), "(for proba = .5)"))
lines(c(0,1), c(0,1))

# Agreements
linkedRecords = which(Delta==1, arr.ind=TRUE)
linkedRecordsFN = which(FalseNegative, arr.ind=TRUE)
linkedRecordsTP = which(TruePositive, arr.ind=TRUE)
for (k in 1:length(pivs))
{
  print(sprintf("PIV %s", pivs[k]))
  agreementD = sum(A[linkedRecords[,1], pivs[k]] == B[linkedRecords[,2], pivs[k]]) / nrow(linkedRecords)
  agreementFN = sum(A[linkedRecordsFN[,1], pivs[k]] == B[linkedRecordsFN[,2], pivs[k]]) / nrow(linkedRecordsFN)
  agreementTP = sum(A[linkedRecordsTP[,1], pivs[k]] == B[linkedRecordsTP[,2], pivs[k]]) / nrow(linkedRecordsTP)
  missingsAD = sum(A[linkedRecords[,1], pivs[k]] == 0) / nrow(linkedRecords)
  missingsBD = sum(B[linkedRecords[,2], pivs[k]] == 0) / nrow(linkedRecords)
  missingsAFN = sum(A[linkedRecordsFN[,1], pivs[k]] == 0) / nrow(linkedRecordsFN)
  missingsBFN = sum(B[linkedRecordsFN[,2], pivs[k]] == 0) / nrow(linkedRecordsFN)
  missingsATP = sum(A[linkedRecordsTP[,1], pivs[k]] == 0) / nrow(linkedRecordsTP)
  missingsBTP = sum(B[linkedRecordsTP[,2], pivs[k]] == 0) / nrow(linkedRecordsTP)
  print(sprintf("FN  agree at %s percent    - FN  missing in A: %s percent    - FN  missing in B: %s percent", round(agreementFN*100,2), round(missingsAFN*100,2), round(missingsBFN*100,2)))
  print(sprintf("D   agree at %s percent    - D   missing in A: %s percent    - D   missing in B: %s percent", round(agreementD*100,2), round(missingsAD*100,2), round(missingsBD*100,2)))
  print(sprintf("TP  agree at %s percent    - TP  missing in A: %s percent    - TP  missing in B: %s percent", round(agreementTP*100,2), round(missingsATP*100,2), round(missingsBTP*100,2)))
}
