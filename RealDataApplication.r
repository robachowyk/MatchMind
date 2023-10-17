################################################################################
#                         NLCTS REAL DATA APPLICATION                         #
################################################################################

library(mvtnorm)
library(numDeriv)
library(Matrix)
library(MASS)
library(plyr)
library(testit)
library(Metrics)
library(haven)
library(tidyr)
library(labelled)
library(progress)

getwd()
dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

A = read.table(file.path("data", "A1982.txt"))
head(A)
summary(A)

B = read.table(file.path("data", "B1994.txt"))
head(B)
summary(B)

regTimeA = A$time
regTimeB = B$time

pivs = c("sex", "dob_yy", "dob_mm", "dob_dd", "reg", "state")
pivs_stable = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)

Linked_id = intersect(A$seq, B$seq)
Nlinks = length(Linked_id)

print("Size of the overlapping set:")
print(Nlinks)
print("Represented as a fraction of the smallest file:")
print(Nlinks / min(nrow(A), nrow(B)))

nvalues = rep(0, length(pivs))
all = rbind( A[,pivs], B[,pivs])
for (k in 1:length(pivs)){
  unique_values = unique(all[,k])
  unique_values = unique_values[unique_values!=0]
  nvalues[[k]] = length(unique_values)
}

nvalues # c(2, 36, 12, 31, 12, 47)

# Add covariates for unstable PIVs survival model
pSameH.varA = list(c(), c(), c(), c(), c("helpertype", "helperpresent"), c("helpertype", "helperpresent"))
pSameH.varB = list(c(), c(), c(), c(), c("helpertype", "helperdo1", "helperdo2"), c("helpertype", "helperdo1", "helperdo2"))
XA = data.frame(A)
XB = data.frame(B)

################################################################################
#                      LAUNCH THE ALGORITHM ON REAL DATA                       #
################################################################################

dataNLTCS = list( pivsA        = A[,pivs], 
                  pivsB        = B[,pivs], 
                  nvalues      = nvalues,
                  pivs_stable  = pivs_stable, 
                  regTimeA     = regTimeA, 
                  regTimeB     = regTimeB,
                  XA           = XA,
                  XB           = XB,
                  pSameH.varA  = pSameH.varA,
                  pSameH.varB  = pSameH.varB)

source("recordlinkage.r")
Rcpp:::sourceCpp("functions.cpp")

fit = stEM( data       = dataNLTCS, 
            nIter      = 250, 
            nBurnin    = 25, 
            MStepIter  = 25, 
            trace      = 1 )

################################################################################
#                            VISUALISE THE RESULTS                             #
################################################################################

par(mfrow=c(1,2))
Delta = matrix(0, nrow=nrow(dataNLTCS$pivsA), ncol=nrow(dataNLTCS$pivsB))
for (lr in 1:Nlinks)
{
  i = which(A$seq==Linked_id[lr])
  j = which(B$seq==Linked_id[lr])
  Delta[i,j] = 1
}
image(x=1:nrow(dataNLTCS$pivsB), y=1:nrow(dataNLTCS$pivsA), z=t(as.matrix(Delta)), xlab="obs. in B", ylab="obs. in A")
title(main = c("True linkage matrix", "underlying the data", sprintf("%s linked record pairs", sum(Delta))), font.main = 4)
image(x=1:nrow(dataNLTCS$pivsB), y=1:nrow(dataNLTCS$pivsA), z=t(as.matrix(fit$Delta)), xlab="obs. in B", ylab="obs. in A")
title(main = c("Estimated linkage matrix fitting the data","taking into account unstability in PIVs", sprintf("%s linked record pairs (proba > .5)", sum(fit$Delta>0.5))), font.main = 4)

par(mfrow=c(1,1))
# gamma
gamma = fit$gamma
plot(gamma, type="l", ylim=c(0,1), xlab = "MC-EM Iterations", ylab = "gamma")
abline(h = Nlinks/min(nrow(A), nrow(B)), col="red")

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
  if(!dataNLTCS$pivs_stable[k]){
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
