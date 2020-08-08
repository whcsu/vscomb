# Based on https://github.com/wwrechard/screening
library(splines)
library(MASS)
library(Boruta)

library(glmnet)

setwd("F:/workingpapers/vscomb/code/results/real")
##

auc<- function(outcome, proba){
  N = length(proba)
  N_pos = sum(outcome)
  df = data.frame(out = outcome, prob = proba)
  df = df[order(-df$prob),]
  df$above = (1:N) - cumsum(df$out)
  return( 1- sum( df$above * df$out ) / (N_pos * (N-N_pos) ) )
}


########################################################################################################################
# Functions to compute the criteria in the simulation.

# Overall selection proportion
PA <- function(activetrue,active)
{
  count <- 0;
  if (setequal(activetrue,intersect(active,activetrue))==TRUE){count <- 1;}
  return(count);
}
# Individual selection proportion
PI <- function(activetrue,active)
{
  q <- length(activetrue);
  count <- rep(0,q);
  for (j in 1:q)
  {
     if (is.element(activetrue[j], active)==TRUE) {count[j] <- 1;}
  }
  return(count);
} 

# 1. M: to compute the minimum model size to ensure the inclusion of all active predictors. 
# 2. mqtl: to compute the 5%, 25%, 50%, 75% and 95% quantiles of the minimum model size out of 1,000 replications.
# 3. Sel.rate: to compute the proportion that every single active predictor is selected 
#    for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.

M<-function(true.v,rank.mtx) {
  # Input
  # true.v   :  the true variables index
  # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
  #             each column corresponds the ranked index in one replication.
  # Output
  # M        :  a vector of the minimum model sizes to ensure the inclusion of all active predictors 
  r<-min(dim(rank.mtx)[2],length(rank.mtx))
  M<-c()
  for (j in 1:r) {M[j]<-max(match(true.v,rank.mtx[,j]))}
  return(M)
}
mqtl<-function(M) {
  # Input    
  # M        :  a vector of the minimum model sizes to ensure the inclusion of all active predictors 
  # Output
  # 5%,25%,50%,75%,95% quantiles of minimum model sizes out of 1000 replications
  quantile(M, probs =c(0.05,0.25,0.5,0.75,0.95))
}

Sel.rate<-function(n,activetrue,rank.mtx) {
  # Input
  # n        :  the sample size
  # c        :  coeficient of cutoffs, for example c=2, cutoff=2[n/log(n)]
  # true.v   :  the true variables index
  # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
  #             each column corresponds the ranked index in one replication.
  # Output
  # rate     :  the proportions that every single active predictor is selected 
  #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
  r<-min(dim(rank.mtx)[2],length(rank.mtx))
  p0<-length(activetrue)
  R<-matrix(0,p0,r)
  rate<-c()
  for (i in 1:p0) {
    for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx[,j]-activetrue[i]))==0) }
    rate[i]<-mean(R[i,])
  }
  return(rate)
}

# Overall selection proportion
PA.rate <- function(activetrue,rank.mtx)
{
  
  r<-min(dim(rank.mtx)[2],length(rank.mtx))
  PA=rep(0,r)
  for (j in 1:r) {
    if (setequal(activetrue,intersect(rank.mtx[,j],activetrue))==TRUE)
        {
      PA[j]=1
    }
    
    }
  return(sum(PA)/r)

  
}

# Proportion Coverage

PC.rate <- function(activetrue,rank.mtx)
{
  
  r<-min(dim(rank.mtx)[2],length(rank.mtx))
  p0<-length(activetrue)
  R<-rep(0,r)

    for (j in 1:r) {
      R[j]<-length(intersect(rank.mtx[,j],activetrue))/length(activetrue)  
      }
    
 
  return (R)
  
  
}


creat.sigma1 <- function(rho,p) {
  Sigma1<-matrix(0,p,p)
  for (i in 1:p){
    for (j in max(1,(i-50)):min(p,(i+50))){
      Sigma1[i,j]<-rho^(abs(i-j))
    }
  }	
  return(Sigma1)
}
#sigma = creat.sigma1(rho,p)
########################################################################################################################

orderedunion<-function(set1,set2,newlen=length(set1)){
  n1=length(set1)  
  n2=length(set2)
  newset=list()
  set2pos=1
  for (set1pos in 1:n1){
    newelm=set1[set1pos]
    
    if  (!is.element(newelm,newset))
    {   
      if (newelm==set2[set2pos]) { 
        newset=append(newset,newelm)
        set2pos=set2pos+1
      }
      else  {
        newset=append(newset,newelm)
        if  (!is.element(set2[set2pos],newset)) 
          newset=append(newset,set2[set2pos])
        set2pos=set2pos+1
      }
    }
    
  }
  for  (i in set2pos:n2) 
    if  (!is.element(set2[i],newset)) newset=append(newset,set2[i])
  newset=unlist(newset)
  
  return (newset[1:newlen])
}


#' An efficient variable screening method
#'
#' This function implements 4 different screening methods (SIS, HOLP, RRCS and Forward regression) for linear models and 3 (excluding RRCS) for generalized linear models.
#' @param x the predictor variables, each row corresponds to an observation. Should be a numeric matrix instead of a data.frame
#' @param y the observation.
#' @param method the screening method to use. Choices are "sis", "holp", "rrcs", "forward". Default to "holp".
#' @param num.select the number of variables to keep after screening. Default to half of the sample size. It will not be used if ebic is set to be TRUE.
#' @param family the model type choices are the same as glmnet. Default to be 'gaussian'.
#' @param ebic Indicate whether the extended BIC should be used to determine the number of variables to keep. If ebic is TRUE, then the algorithm will use ebic to terminate the screening procedure and num.select will be ignored.
#' @param ebic.gamma tunning parameter for ebic (between 0 and 1). Gamma = 0 corresponds to the usual BIC. default to be 1.
#' @return a list of two variables "screen" and "method". "screen" contains the index of the selected variables and "method" indicates the method of the screening.
#'
#' @examples There are one unit test function and two integrated test functions. Two integrated function test on linear model and logistic model. User specify the sample size, dimension and the true indexes. The two function generate simulate data and coefficients and print the screening results for all methords.
#'
#' linearModelTest(n = 50, p = 100, beta.not.null = c(1, 2, 3), num.select = 20)
#' logisticTest(n = 50, p = 100, beta.not.null = c(1, 2, 3), nums.select = 20)
#'
#' @references
#' Fan, Jianqing, and Jinchi Lv. "Sure independence screening for ultrahigh dimensional feature space." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70.5 (2008): 849-911.
#' Wang, Xiangyu, and Chenlei Leng. "High-dimensional ordinary least-squares projection for screening variables." arXiv preprint arXiv:1506.01782 (2015).
#' Li, Gaorong, et al. "Robust rank correlation based screening." The Annals of Statistics 40.3 (2012): 1846-1877.
#' Wang, Hansheng. "Forward regression for ultra-high dimensional variable screening." Journal of the American Statistical Association 104.488 (2009): 1512-1524.
#'
screening <- function(x, y, method = 'holp', num.select = floor(dim(x)[1]/log(dim(x)[1])), family = 'gaussian', ebic = FALSE, ebic.gamma = 1) {
  
  # standardize
  x = as.matrix(x)
  X = scale(x)
  if (family == 'gaussian'){
    Y = y - mean(y)
  }
  else {
    Y = y
  }
  if (is.null(dim(X))) {
    p = 1
    n = length(X)
  } else {
    n = dim(X)[1]
    p = dim(X)[2]
  }
  
  # if p is smaller than the required number of variables, return all
  if (p == 1 || (p < num.select && !ebic)) {
    selectedVariable = 1 : p
  }
  else {
    # for the linear case, it is easy to compute everything.
    if (family == 'gaussian') {
      if (method == 'holp') {
        OLS = t(X) %*% solve(X %*% t(X) + diag(n) * 1, Y)
        ranking = sort(abs(OLS), index.return = TRUE, decreasing = TRUE)
        if (ebic) {
          result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
        }
        else {
          result = ranking$ix[1 : num.select]
        }
      }
      else if (method == 'sis') {
        OLS = t(X) %*% Y
        ranking = sort(abs(OLS), index.return = TRUE, decreasing = TRUE)
        if (ebic) {
          result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
        }
        else {
          result = ranking$ix[1 : num.select]
        }
      }
      else if (method == 'rrcs') {
        OLS = .rankScreening(X, Y)
        ranking = sort(abs(OLS), index.return = TRUE, decreasing = TRUE)
        if (ebic) {
          result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
        }
        else {
          result = ranking$ix[1 : num.select]
        }
      }
      else if (method == 'forward') {
        if (ebic) {
          ranking = .forwardRegression(X, Y, n - 1, family, ebic.gamma)
          bestModel = which.min(ranking$bic)
          result = ranking$select[1 : bestModel]
        }
        else {
          result = .forwardRegression(X, Y, num.select, family, ebic.gamma)$select
        }
      }
      else {
        stop('The method is unknown and not supported/')
      }
    }
    else {
      if (method == 'holp') {
        require(glmnet)
        model = glmnet(x = X, y = Y, family = family, alpha = 0, lambda = 1, intercept = FALSE)
        coefs = coef(model)[-1]
        ranking = sort(abs(coefs), index.return = TRUE, decreasing = TRUE)
        if (ebic) {
          result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
        }
        else {
          result = ranking$ix[1 : num.select]
        }
      }
      else if (method == 'sis') {
        residuals = rep(0, p)
        for (i in 1 : p) {
          model = glm(Y ~ X[, i] - 1, family = family)
          residuals[i] = model$deviance
        }
        ranking = sort(residuals, index.return = TRUE)
        if (ebic) {
          results = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
        }
        else {
          result = ranking$ix[1 : num.select]
        }
      }
      else if (method == 'forward') {
        if (ebic) {
          ranking = .forwardRegression(X, Y, n - 1, family, ebic.gamma)
          bestModel = which.min(ranking$bic)
          result = ranking$select[1 : bestModel]
        }
        else {
          result = .forwardRegression(X, Y, num.select, family, ebic.gamma)$select
        }
      }
      else if (method == 'rrcs') {
        stop('rrcs is not supported for GLM')
      }
      else {
        stop('The method is unknown and not supported.')
      }
    }
  }
  return (list(screen = result, method = method))
}

#############################################################################

.rankScreening <- function(X, Y) {
  n = dim(X)[1]
  p = dim(X)[2]
  w = rep(0,p)
  Ynew = sort(Y, index.return = T)
  for(j in 1 : n)
    for(k in j : n)
      w = w + (X[Ynew$ix[k], ] > X[Ynew$ix[j], ])
  w = w / n / (n-1)
  w = w - 1 / 4
  return (abs(w))
}

#############################################################################

.forwardRegression <- function(X, Y, num.select, family, ebic.gamma = 1) {
  
  # get the dimension
  n = dim(X)[1]
  p = dim(X)[2]
  
  # store the used variables, including good variable and bad variable
  usedVariables = NULL
  
  # store the selected variables
  selectedVariables = NULL
  
  # store the residual sum of squares
  rss = rep(0, p)
  
  # store the bic values
  bic = rep(Inf, num.select)
  
  iteration = 0
  while (iteration < num.select) {
    
    # to compute for each variable the deviance of the model if they were added
    for (i in setdiff(1 : p, usedVariables)) {
      activeVariables = c(selectedVariables, i)
      #-1 denotes without intercept
      model = try(glm(Y ~ X[, activeVariables] - 1, family = family))
      if (inherits(model, "try-error")) {
        rss[i] = Inf
        usedVariables = c(usedVariables, i)
      } else {
        rss[i] = model$deviance
      }
    }
    if (min(rss) == Inf) {
      break
    }
    
    # select the variabel that gives the smallest deviance to the model
    toAdd = which.min(rss)
    selectedVariables = c(selectedVariables, toAdd)
    usedVariables = c(usedVariables, toAdd)
    
    # record the corresponding bic value
    iteration = iteration + 1
    bic[iteration] = .ebic(rss[toAdd], p, n, iteration, ebic.gamma)
    rss[toAdd] = Inf
  }
  
  return (list(select = selectedVariables, bic = bic))
}

##############################################################################

.ebic <- function(deviance, model.size, sample.size, num.select, ebic.gamma) {
  return (deviance + num.select * (log(sample.size) + 2 * ebic.gamma * log(model.size)))
}

##############################################################################

.ebicRanking <- function(X, Y, sortedVariables, family, ebic.gamma) {
  # get the dimension
  n = dim(X)[1]
  p = dim(X)[2]
  
  # store the currently selected variables
  selectedVariables = NULL
  
  # store the bic values
  bic = rep(Inf, n - 1)
  
  iteration = 0
  while (iteration < n - 1) {
    
    iteration = iteration + 1
    
    # to compute for each variable the deviance of the model if they were added
    i = sortedVariables[iteration]
    selectedVariables = c(selectedVariables, i)
    model = try(glm(Y ~ X[, selectedVariables] - 1, family = family))
    if (inherits(model, "try-error")) {
      rss = Inf
    } else {
      rss = model$deviance
    }
    
    if (rss == Inf) {
      break
    }
    
    # record the corresponding bic value
    bic[iteration] = .ebic(rss, p, n, iteration, ebic.gamma)
  }
  
  bestModel = which.min(bic)
  return (list(select = sortedVariables[1 : bestModel], bic = bic))
}

set.seed(123)

# 5 folds 10 times, total running time 50 times

Rn=10
totalfold=2

for (i in 19:19){

datafun=paste0("mydata",i)
mydat=eval(parse(text=datafun))()

mydata=mydat$curdata
rii=mydat$currii
datasetname=mydat$curdatasetname


x=mydata[,-c(rii)]
y=mydata[,c(rii)]


num_select=floor(dim(x)[1]/log(dim(x)[1]))
#if(num_select%%2!=0) num_select=num_select+1

ss1=screening(x,y,num.select=num_select,method="sis")
ss2=screening(x,y,num.select=num_select,method="rrcs")
ss3=screening(x,y,num.select=num_select,method="holp")
rfs=Boruta(x,y); 
ssall=order(attStats(rfs)$medianImp,decreasing = T)
RF=ssall[1:num_select]
#ss4_half=ssall[1:(num_select/2)]

SIS=ss1$screen[1:(num_select)]
RRCS=ss2$screen[1:(num_select)]
HOLP=ss3$screen[1:(num_select)]

R_S=orderedunion(ss1$screen[1:(num_select)],RF)
R_R=orderedunion(ss2$screen[1:(num_select)],RF)
R_H=orderedunion(ss3$screen[1:(num_select)],RF)

lassoauc=c(rep(0,Rn))
sisLasso=c(rep(0,Rn))
RRCSLasso=c(rep(0,Rn))
HOLPLasso=c(rep(0,Rn))
RFLasso=c(rep(0,Rn))
R_SLasso=c(rep(0,Rn))
R_RLasso=c(rep(0,Rn))
R_HLasso=c(rep(0,Rn))
testresult <-data.frame(row.names=1:Rn)


n=dim(mydata)[1]


#Randomly shuffle the data
for (ii in seq(1, Rn, by=totalfold)){
  print(ii)
  
  mydata<-mydata[sample(nrow(mydata)),]
  
  #Create 5 equally size folds
  folds <- cut(seq(1,nrow(mydata)),breaks=totalfold,labels=FALSE)
  #Perform 2 fold cross validation
  
  for(k in 1:totalfold){
    #Segement your data by fold using the which() function
    testIndexes <- which(folds==k,arr.ind=TRUE)
    tesety <- mydata[testIndexes, rii]
    trsety <- mydata[-testIndexes, rii]
    
    tesetx <- mydata[testIndexes, -c(rii)]
    trsetx <- mydata[-testIndexes, -c(rii)]
    
    glmfit=cv.glmnet(x=as.matrix(trsetx),y=trsety,alpha=1,family="binomial")
    glmpre=predict(glmfit,as.matrix(tesetx), type="response",s="lambda.min")
    lassoauc[ii+k-1]=auc(tesety,as.vector(glmpre))
    
    cat("\nLasso\n",lassoauc[ii+k-1])
    
    
    sis_glmfit=cv.glmnet(x=as.matrix(trsetx[,SIS]),y=trsety,alpha=1,family="binomial")
    
    sis_glmpre=predict(sis_glmfit,as.matrix(tesetx[,SIS]), type="response",s="lambda.min")
    
    sisLasso[ii+k-1]=auc(tesety,as.vector(sis_glmpre))
    cat("\nsisLasso\n",sisLasso[ii+k-1])
    
    
    RRCS_glmfit=cv.glmnet(x=as.matrix(trsetx[,RRCS]),y=trsety,alpha=1,family="binomial")
    
    RRCS_glmpre=predict(RRCS_glmfit,as.matrix(tesetx[,RRCS]), type="response",s="lambda.min")
    
    RRCSLasso[ii+k-1]=auc(tesety,as.vector(RRCS_glmpre))
    
    cat("\nRRCSLasso\n",RRCSLasso[ii+k-1])
    
    
    HOLP_glmfit=cv.glmnet(x=as.matrix(trsetx[,HOLP]),y=trsety,alpha=1,family="binomial")
    
    HOLP_glmpre=predict(HOLP_glmfit,as.matrix(tesetx[,HOLP]), type="response",s="lambda.min")
    
    HOLPLasso[ii+k-1]=auc(tesety,as.vector(HOLP_glmpre))
    
    cat("\nHOLPLasso\n",HOLPLasso[ii+k-1])
    
    
    RF_glmfit=cv.glmnet(x=as.matrix(trsetx[,RF]),y=trsety,alpha=1,family="binomial")
    
    RF_glmpre=predict(RF_glmfit,as.matrix(tesetx[,RF]), type="response",s="lambda.min")
    
    RFLasso[ii+k-1]= auc(tesety,as.vector(RF_glmpre))
    
    cat("\nRFLasso\n",RFLasso[ii+k-1])
    
    
    
    R_S_glmfit=cv.glmnet(x=as.matrix(trsetx[,R_S]),y=trsety,alpha=1,family="binomial")
    
    R_S_glmpre=predict(R_S_glmfit,as.matrix(tesetx[,R_S]), type="response",s="lambda.min")
    
    R_SLasso[ii+k-1]= auc(tesety,as.vector(R_S_glmpre))
    
    cat("\nR_SLasso\n",R_SLasso[ii+k-1])
    
    
    R_R_glmfit=cv.glmnet(x=as.matrix(trsetx[,R_R]),y=trsety,alpha=1,family="binomial")
    
    R_R_glmpre=predict(R_R_glmfit,as.matrix(tesetx[,R_R]), type="response",s="lambda.min")
    
    R_RLasso[ii+k-1]=auc(tesety,as.vector(R_R_glmpre))
    
    cat("\nR_RLasso\n",R_RLasso[ii+k-1])
    
    
    R_H_glmfit=cv.glmnet(x=as.matrix(trsetx[,R_H]),y=trsety,alpha=1,family="binomial")
    
    R_H_glmpre=predict(R_H_glmfit,as.matrix(tesetx[,R_H]), type="response",s="lambda.min")
    
    R_HLasso[ii+k-1]=auc(tesety,as.vector(R_H_glmpre))
    
    cat("\nR_HLasso\n",R_HLasso[ii+k-1])
    
  }
  
}
#setwd("")
testresult=cbind(testresult,data.frame(lassoauc,sisLasso,RRCSLasso,HOLPLasso,RFLasso,R_SLasso,R_RLasso,R_HLasso))
csvfn=paste0("5_2fold_",datasetname,".csv")
write.csv(testresult,csvfn)
boxplot(testresult)
}

