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
library(beepr)

getwd()
dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

################################################################################

A = read.table(file.path("data", "A1982.txt"))
summary(A)

B = read.table(file.path("data", "B1994.txt"))
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

nvalues = c(2, 36, 12, 31, 12, 47)

# ADD COVARIATES FOR UNSTABLE PIVs SURVIVAL MODELs
pSameH.varA = list(c(), c(), c(), c(), c("helpertype", "helperpresent"), c("helpertype", "helperpresent"))
pSameH.varB = list(c(), c(), c(), c(), c("helpertype", "helperdo1", "helperdo2"), c("helpertype", "helperdo1", "helperdo2"))
XA = data.frame(A)
XB = data.frame(B)

dataNLTCS= list(pivsA = A[,pivs], 
                pivsB = B[,pivs], 
                nvalues = nvalues,
                pivs_stable = pivs_stable, 
                regTimeA = regTimeA, 
                regTimeB = regTimeB,
                XA = XA,
                XB = XB,
                pSameH.varA = pSameH.varA,
                pSameH.varB = pSameH.varB)

Delta = matrix(0, nrow=nrow(dataNLTCS$pivsA), ncol=nrow(dataNLTCS$pivsB))
for (i in 1:nrow(dataNLTCS$XA))
{
  idA = A$seq[i]
  if(idA %in% Linked_id)
  {
    j = which(B$seq==idA)
    Delta[i,j] = 1
  }
}
par(mfrow=c(1,2))
image(x=1:nrow(dataNLTCS$pivsB), y=1:nrow(dataNLTCS$pivsA), z=t(as.matrix(Delta)), xlab="obs. in B", ylab="obs. in A") # gray.colors(10)
title(main = c("True linkage matrix", "underlying the data"), font.main = 4)

pdf(file = "LinkageTruth.pdf")
image(x=1:nrow(dataNLTCS$pivsB), y=1:nrow(dataNLTCS$pivsA), z=t(as.matrix(Delta)), xlab="obs. in B", ylab="obs. in A") # gray.colors(10)
title(main = c("True linkage matrix", "underlying the data"), font.main = 4)
dev.off()

source("recordlinkage_new.r")
Rcpp:::sourceCpp("functions_new.cpp")

# A SONG WILL START AT THE END SO YOU CAN DO SOMETHING ELSE WHILE WAITING
start = Sys.time()
print("Launch RL algorithm:")
fit = stEM(data = dataNLTCS, 
           nIter = 125, 
           nBurnin = 20, 
           MStepIter = 25, 
           trace = 1)
end = Sys.time()
print(end-start)
browseURL('https://www.youtube.com/watch?v=NTa6Xbzfq1U')

# PLOT FOR GAMMA
par(mfrow=c(2,3))
gamma = fit$gamma
plot(gamma, type="l", ylim=c(0,1), xlab = "MC-EM Iterations", ylab = "gamma")
abline(h = Nlinks / min(nrow(A), nrow(B)), col="red")

# PLOT FOR ETA
par(mfrow=c(2,3))
eta = fit$eta
for(k in 1:length(eta))
{
  vec = c("eta 1")
  plot(eta[[k]][,1], ylim=c(0,1), type="l", col=1, xlab = "MC-EM Iterations", ylab = sprintf("PIV %s", k))
  for(j in 2:ncol(eta[[k]]))
  {
    lines(eta[[k]][,j], ylim=c(0,1), col=j)
    vec <- append(vec, sprintf("eta %s", j))
  }
  legend("topright", legend=vec, col=seq_len(ncol(eta[[k]])), lty=1:4, cex=0.8)
} 

# PLOT FOR OMEGA
par(mfrow=c(2,3))
omega = fit$omega
for(k in 1:length(omega))
{
  if(!dataNLTCS$pivs_stable[k]){
    vec = c("time difference")
    plot(omega[[k]][,1], ylim=c(min(omega[[k]])-0.5, max(omega[[k]])+0.5), type="l", col=1, xlab = "MC-EM Iterations", ylab = sprintf("Survival model coef. for unstable PIV %s", k))
    if(ncol(omega[[k]])>2){
      vec <- append(vec, dataNLTCS$pSameH.varA[[k]])
      vec <- append(vec, dataNLTCS$pSameH.varB[[k]])
      for(c in 2:ncol(omega[[k]]))
      {
        lines(omega[[k]][,c], col=c)
      }
    }
    legend("topright", legend=vec, col=c(1:ncol(omega[[k]])), lty=1:4, cex=0.8)
  }
}

# PLOT FOR PHI
par(mfrow=c(2,3))
phi = fit$phi
for(k in 1:length(phi))
{
  plot(phi[[k]][,1], ylim=c(0,1), type="l", col=1, xlab = "MC-EM Iterations", ylab = sprintf("PIV %s", k))
  lines(phi[[k]][,2], col=2)
  lines(phi[[k]][,3], col=3)
  legend("topright", legend=c("agree", "miss A", "miss B"), col=c(1,2,3), lty=1:4, cex=0.8)
}

# CONFUSION MATRIX AND ROC-AUC
# if negative: fit_unstable$Delta>0.5 < Delta: should be a match but is not: FN
# if positive: fit_unstable$Delta>0.5 > Delta: should not be a match but is: FP
# if neutral: fit_unstable$Delta>0.5 = Delta: if Delta is 1: TP if Delta is 0: TN
DeltaSize = nrow(A)*nrow(B)
sprintf("matrix size: %s", DeltaSize)
sprintf("nbr of links: %s", Nlinks)

TP = sum( ((fit$Delta>0.5) - Delta == 0) & (Delta == 1) )
TN = sum( ((fit$Delta>0.5) - Delta == 0) & (Delta == 0) )
FP = sum((fit$Delta>0.5) - Delta > 0)
FN = sum((fit$Delta>0.5) - Delta < 0)
FPR = FP / (FP + TN)
TPR = TP / (TP + FN)
sprintf("TP: %s", TP)
sprintf("TN: %s", TN)
sprintf("FP: %s", FP)
sprintf("FN: %s", FN)

par(mfrow=c(2,2))
FPR_list = c()
TPR_list = c()
for(t in seq(0, 1.1, by=0.1)){
  TP = sum( ((fit$Delta>t) - Delta == 0) & (Delta == 1) )
  TN = sum( ((fit$Delta>t) - Delta == 0) & (Delta == 0) )
  FP = sum((fit$Delta>t) - Delta > 0)
  FN = sum((fit$Delta>t) - Delta < 0)
  FPR_list = append(FPR_list, FP / (FP + TN))
  TPR_list = append(TPR_list, TP / (TP + FN))
}
plot(FPR_list, TPR_list, xlim=c(0,1), ylim=c(0,1))
title(sprintf("AUC: %s", auc(Delta, fit$Delta)))
lines(c(0,1), c(0,1))
plot(FPR_list, TPR_list, xlim=c(0,0.002), ylim=c(0,1))
lines(c(0,1), c(0,1))
