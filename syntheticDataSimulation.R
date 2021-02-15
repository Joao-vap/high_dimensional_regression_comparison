library(MASS)                       
library(glmnet)
#library(parcor)
library(ncvreg)
library(monomvn)
library(flare)                  
library(c060)
library(Matrix)
library(ROCR)

################################################
# data generation function
################################################

data_generation <- function(n, p, s0, SNR, 
                            corDesign = c("Independence", "Pairwise", "Toeplitz"), 
                            pB = 0, rho = 0, s0B = 0,
                            n_test = 500) {
  
  corDesign <- match.arg(corDesign)
  n_training=n
  n=n_training+n_test  
  
  
  if(corDesign=="Independence") {
    X=matrix(rnorm(n*p),n,p)   
    gamma=rep(0,p); gamma[sample(1:p,size=s0,replace=FALSE)]=1
    beta0=rep(0,p); beta0[which(gamma==1)]=beta_value
  }
  if(corDesign=="Pairwise") {
    Block_rho=rho
    n_block=p/pB
    block=function() {
      signal_idx=sample(1:pB,size=s0B,replace=FALSE)   ######Where you put the signals in each block#####
      e0=rep(1,pB)          #####A vector of 1's####
      Sigma=Block_rho*(e0%*%t(e0)-diag(pB))+diag(pB)
      result=mvrnorm(n=n,mu=rep(0,pB),Sigma=Sigma,tol=1e-12)    #####mvrnorm has realization in each ROW#######
      colnames(result)=rep(0,pB)
      colnames(result)[signal_idx]=1
      result
    }      
    
    X=block(); if(n_block>=2) {for(b in 1:(n_block-1)){X=cbind(X,block())}}
    gamma=as.numeric(colnames(X)); beta0=rep(0,p); beta0[which(gamma==1)]=beta_value
    
  }     
  
  if(corDesign=="Toeplitz") {
    rho_base=0.95
    n_block=p/pB
    cor_between_two_signals=rho
    no_signal_per_block=s0B
    
    block=function() {
      Sigma=matrix(0,pB,pB)
      for(i in 1:ncol(Sigma))                 
      {
        for(j in 1:i)                            
        {Sigma[j,i]=rho_base^(abs(i-j))}
      }
      Sigma[lower.tri(Sigma)] = t(Sigma)[lower.tri(Sigma)]     
      result=mvrnorm(n=n,mu=rep(0,pB),Sigma=Sigma,tol=1e-12)    
      colnames(result)=rep(0,pB) ################################# Block_size?? Maybe pB? ##############################
      if(no_signal_per_block!=2) {stop("Should be 2 signals per block for Toeplitz Design!")}
      dis=round(log(cor_between_two_signals,base=rho_base))
      if(1+dis>=pB) {stop("The block size is too small for Toeplitz Design!")}
      signal_idx=c(round(pB/2),round(pB/2)+dis)      
      colnames(result)[signal_idx]=1  
      result
    }
    X=block(); for(b in 1:(n_block-1)){X=cbind(X,block())}
    gamma=as.numeric(colnames(X)); beta0=rep(0,p); beta0[which(gamma==1)]=beta_value
  }     
  
  
  non_zero_idx=which(beta0!=0); if(length(non_zero_idx)>s0) {beta0[-non_zero_idx[1:s0]]=0} 
  #################for correlation designs, if not all blocks should have signals, those blocks besides first s0 signals are empty##########
  
  
  Xr=X; X_training=Xr[1:n_training,];X_test=Xr[(n_training+1):n,]
  X_training=scale(X_training); X_test=scale(X_test)
  sigmae=as.numeric(sqrt(t(beta0)%*%t(X_training)%*%X_training%*%beta0/(n_training*SNR^2)))
  Xr=rbind(X_training,X_test)              
  Y0=Xr%*%beta0
  
  Y0_training=Y0[1:n_training]
  e_training=rnorm(n_training,0,sigmae)
  Y_training=Y0_training+e_training
  Y_training=Y_training-mean(Y_training)
  Y0_test=Y0[(n_training+1):n]
  e_test=rnorm(n_test,0,sigmae)
  Y_test=Y0_test+e_test
  Y_test=Y_test-mean(Y_test)
  
  
  if(length(e_training)!=length(Y0_training) | length(e_test)!=length(Y0_test)) {stop("WRONG!")}
  
  
  result=list()
  result$X=X_training; result$Y=Y_training; result$X_test=X_test; result$Y_test=Y_test; result$beta0=beta0; result$sigma=sigma
  result
}                    

################################################
# load scenarios
################################################

scenarios <- read.table("syntheticDataScenarios.txt", header = TRUE,  stringsAsFactors = FALSE)

################################################
# run simulation for the selected scenario
################################################

alpha1=0.3; alpha2=0.6; beta_value=3
method_name=c("lasso","lenet","henet","ridge","dantzig","scad","stability")
for(name in method_name) {assign(sprintf("%s_list",name),list())}

G.scenario <- c(); G.method <- c(); G.r <- c(); G.n <- c(); G.p <- c(); G.s0 <- c()
G.SNR <- c(); G.corDesign <- c(); G.rho <- c();  G.pB <- c(); G.s0B <- c()

G.mean_pauc <- c(NA); G.sd_pauc <- c(NA)
G.mean_rmse <- c(NA); G.sd_rmse <- c(NA)
G.mean_tpr <- c(NA); G.sd_tpr <- c(NA)
G.mean_ppv <- c(NA); G.sd_ppv <- c(NA)
G.mean_mcc <- c(NA); G.sd_mcc <- c(NA)

niter=64

for(cen in 1:1899)
  {
  
  ################################################
  # select a scenario
  ################################################
  
  list2env(scenarios[cen, ], .GlobalEnv)
  
  for(iter in 1:64) {
    DATA=data_generation(n, p, s0, SNR, corDesign, pB, rho, s0B)
    X=DATA$X;Y=DATA$Y;X_test=DATA$X_test;Y_test=DATA$Y_test;
    beta0=DATA$beta0; sigma=DATA$sigma
    
    nfolds=5      
    foldid=sample(rep(1:nfolds,n/nfolds),n)
    xlim=c(0,50/(p-s0))
    
    
    #############################################################Lasso###############################################
    lasso=list()
    cv=cv.glmnet(x=X,y=Y,alpha=1,foldid=foldid,nlambda=100)
    ptm=proc.time()
    lasso_run=glmnet(x=X,y=Y,alpha=1,lambda=cv$lambda)
    run_time=proc.time()-ptm
    solution_path=Matrix(lasso_run$beta,sparse=TRUE)
    colnames(solution_path)=cv$lambda
    lasso$run.time=run_time
    lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
    beta.cv=solution_path[,lambda.min_idx]
    tp=sum(beta.cv!=0 & beta0!=0)
    fp=sum(beta.cv!=0 & beta0==0)
    tn=sum(beta.cv==0 & beta0==0)
    fn=sum(beta.cv==0 & beta0!=0)
    if((tp+fp)==0) {mcc=0;ppv=0}
    if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
    tpr=tp/(tp+fn)
    lasso$mcc=mcc
    lasso$ppv=ppv
    lasso$tpr=tpr
    
    
    rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
    lasso$rmse=rmse
    score=rep(0,p)
    for(var_idx in 1:p)
    {
      if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}  ######Never entered solution path#####
      if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}   ######Always in solution path, save the largest lambda####
      
      if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
      {
        idx=max(which(solution_path[var_idx,]==0))
        if(idx!=length(cv$lambda)) {idx=idx+1; score[var_idx]=cv$lambda[idx]}
      }
    }
    pred=prediction(abs(score),ifelse(beta0!=0,1,0))
    auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
    pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
    lasso$pauc=pauc
    lasso_list[[iter]]=lasso
    
    
    ###############################################################Lenet#############################################
    lenet=list()       
    cv=cv.glmnet(x=X,y=Y,alpha=alpha2,foldid=foldid,nlambda=100)
    ptm=proc.time()
    lenet_run=glmnet(x=X,y=Y,alpha=alpha2,lambda=cv$lambda)
    run_time=proc.time()-ptm
    solution_path=Matrix(lenet_run$beta,sparse=TRUE)
    colnames(solution_path)=cv$lambda
    lenet$run.time=run_time
    lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
    beta.cv=solution_path[,lambda.min_idx]
    
    tp=sum(beta.cv!=0 & beta0!=0)
    fp=sum(beta.cv!=0 & beta0==0)
    tn=sum(beta.cv==0 & beta0==0)
    fn=sum(beta.cv==0 & beta0!=0)
    if((tp+fp)==0) {mcc=0;ppv=0}
    if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
    tpr=tp/(tp+fn)
    lenet$mcc=mcc
    lenet$ppv=ppv
    lenet$tpr=tpr
    rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
    lenet$rmse=rmse
    
    score=rep(0,p)
    for(var_idx in 1:p)
    {
      if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}   
      if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}
      if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
      {
        idx=max(which(solution_path[var_idx,]==0))
        if(idx!=length(cv$lambda)) {idx=idx+1;score[var_idx]=cv$lambda[idx]}
      }
    }     
    pred=prediction(abs(score),ifelse(beta0!=0,1,0))
    auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
    pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
    lenet$pauc=pauc
    lenet_list[[iter]]=lenet
    
    
    ##############################################################Henet##################
    henet=list()
    cv=cv.glmnet(x=X,y=Y,alpha=alpha1,foldid=foldid,nlambda=100)
    ptm=proc.time()
    henet_run=glmnet(x=X,y=Y,alpha=alpha1,lambda=cv$lambda)
    run_time=proc.time()-ptm
    solution_path=Matrix(henet_run$beta,sparse=TRUE)
    colnames(solution_path)=cv$lambda
    henet$run.time=run_time
    
    lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
    beta.cv=solution_path[,lambda.min_idx]
    tp=sum(beta.cv!=0 & beta0!=0)
    fp=sum(beta.cv!=0 & beta0==0)
    tn=sum(beta.cv==0 & beta0==0)
    fn=sum(beta.cv==0 & beta0!=0)
    if(tp+fp+tn+fn!=p) {warning("wrong tp/fp/tn/fn!")}
    if((tp+fp)==0) {mcc=0;ppv=0}
    if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
    tpr=tp/(tp+fn)
    henet$mcc=mcc
    henet$ppv=ppv
    henet$tpr=tpr
    rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
    henet$rmse=rmse
    
    
    score=rep(0,p)
    for(var_idx in 1:p)
    {
      if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}   
      if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}
      if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
      {
        idx=max(which(solution_path[var_idx,]==0))
        if(idx!=length(cv$lambda)) {idx=idx+1; score[var_idx]=cv$lambda[idx]}
      }
    }    
    pred=prediction(abs(score),ifelse(beta0!=0,1,0))
    auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
    pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
    henet$pauc=pauc
    henet_list[[iter]]=henet
    
    
    ####################################################################Ridge Estimator#########################################################
    ridge=list()
    cv=cv.glmnet(x=X,y=Y,alpha=0,foldid=foldid,nlambda=100)
    ptm=proc.time()
    ridge_run=glmnet(x=X,y=Y,alpha=0,lambda=cv$lambda)
    run_time=proc.time()-ptm
    
    solution_path=Matrix(ridge_run$beta,sparse=TRUE)
    colnames(solution_path)=cv$lambda
    ridge$run.time=run_time
    
    Beta=predict(ridge_run,type="coef",s=cv$lambda.min)
    beta.cv=Beta[-1]
    rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
    ridge$rmse=rmse
    
    pred=prediction(abs(beta.cv),ifelse(beta0!=0,1,0))
    auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
    pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
    ridge$pauc=pauc
    ridge_list[[iter]]=ridge
    
     
    ####################################################################SCAD####################################################
    scad=list()
    cv=cv.ncvreg(X=X,y=Y,penalty="SCAD",nfolds=nfolds,max.iter=5000,nlambda=100)
    ptm=proc.time()
    scad_run=ncvreg(X=X,y=Y,family="gaussian",penalty="SCAD",max.iter=5000,lambda=cv$lambda)
    run_time=proc.time()-ptm
    
    
    solution_path=Matrix(scad_run$beta[-1,],sparse=TRUE)
    colnames(solution_path)=cv$lambda
    if(nrow(solution_path)!=p) {stop("Something wrong with solution path!")}
    scad$run.time=run_time
    lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
    beta.cv=solution_path[,lambda.min_idx]
    
    tp=sum(beta.cv!=0 & beta0!=0)
    fp=sum(beta.cv!=0 & beta0==0)
    tn=sum(beta.cv==0 & beta0==0)
    fn=sum(beta.cv==0 & beta0!=0)
    if(tp+fp+tn+fn!=p) {warning("wrong tp/fp/tn/fn!")}
    if((tp+fp)==0) {mcc=0;ppv=0}
    if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
    tpr=tp/(tp+fn)
    scad$mcc=mcc
    scad$ppv=ppv
    scad$tpr=tpr
    rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
    scad$rmse=rmse
    
    score=rep(0,p)
    for(var_idx in 1:p)
    {
      if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}   
      if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}
      
      if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
      {
        idx=max(which(solution_path[var_idx,]==0))
        if(idx!=length(cv$lambda)) {idx=idx+1; score[var_idx]=cv$lambda[idx]}
      }
    }     
    pred=prediction(abs(score),ifelse(beta0!=0,1,0))
    auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
    pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
    scad$pauc=pauc
    scad_list[[iter]]=scad
    
     
    ######################################################Dantzig Selector#####################################
    
    # dantzig=list()
    # nlambda=100    ######number of lambda we try for finding the optimal one#####
    # ptm=proc.time()
    # dantzig_run=slim(X=X,Y=Y,method="dantzig",nlambda=nlambda,verbose=FALSE,lambda.min.ratio=0.01)
    # run_time=proc.time()-ptm
    # 
    # solution_path=Matrix(dantzig_run$beta,sparse=TRUE)
    # colnames(solution_path)=dantzig_run$lambda
    # 
    # #############k-fold cross-validation for each value of lambda_seq######
    # lambda_seq=dantzig_run$lambda
    # MSE_SEQ=c()   #####mean cross-validated MSE for each possible regularization level (lambda)######
    # MSE_SD_SEQ=c() #####Corresponding standard errors##########
    # k=nfolds  #####number of folds for cross-validation#####
    # samples=sample(1:n,n,replace=FALSE)
    # 
    # max=round(n/k)
    # x=seq_along(samples)
    # groups=split(samples,ceiling(x/max))   #####k groups(equal sizes) samples#####
    # for(i in 1:length(groups)){assign(sprintf("idx_%s",i),groups[[i]])}
    # 
    # for(j in 1:length(lambda_seq))
    # {
    #   lambda=lambda_seq[j]
    #   mse=c()
    #   for(k in 1:length(groups))
    #   {
    #     idx_test=get(sprintf("idx_%s",k))
    #     idx_train=setdiff(1:n,idx_test)
    #     X_test.part=X[idx_test,]
    #     Y_test.part=Y[idx_test]
    #     X_train.part=X[idx_train,]
    #     Y_train.part=Y[idx_train]
    # 
    #     dantzig_train=slim(X=X_train.part,Y=Y_train.part,method="dantzig",lambda=lambda,verbose=FALSE)
    #     beta_train=dantzig_train$beta
    #     Y_test_hat=X_test.part%*%beta_train+rep(dantzig_train$intercept,nrow(X_test.part))
    #     mse=c(mse,mean((Y_test.part-Y_test_hat)^2))
    #   }
    #   MSE_SEQ[j]=mean(mse)
    #   MSE_SD_SEQ[j]=sd(mse)
    # }
    # 
    # 
    # lambda.min=lambda_seq[match(min(MSE_SEQ),MSE_SEQ)]
    # dantzig$run.time=run_time
    # 
    # lambda.min_idx=which.min(abs(lambda_seq-lambda.min))
    # beta.cv=solution_path[,lambda.min_idx]
    # dantzig$beta.cv=beta.cv
    # 
    # tp=sum(beta.cv!=0 & beta0!=0)
    # fp=sum(beta.cv!=0 & beta0==0)
    # tn=sum(beta.cv==0 & beta0==0)
    # fn=sum(beta.cv==0 & beta0!=0)
    # if((tp+fp)==0) {mcc=0;ppv=0}
    # if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
    # tpr=tp/(tp+fn)
    # dantzig$mcc=mcc
    # dantzig$ppv=ppv
    # dantzig$tpr=tpr
    # rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
    # dantzig$rmse=rmse
    # 
    # score=rep(0,p)
    # for(var_idx in 1:p)
    # {
    #   if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}
    #   if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=lambda_seq[1]}
    #   if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
    #   {
    #     idx=max(which(solution_path[var_idx,]==0))
    #     if(idx!=length(lambda_seq)) {idx=idx+1; score[var_idx]=dantzig_run$lambda[idx]}
    #   }
    # }
    # pred=prediction(abs(score),ifelse(beta0!=0,1,0))
    # auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
    # pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
    # dantzig$pauc=pauc
    # dantzig_list[[iter]]=dantzig
    
    
    ######################################################Stability Selection#################################
    stability=list()
    ptm=proc.time()
    stability_run=stabpath(y=Y,x=X,size=0.632,steps=500,weakness=1,nlambda=100)
    run_time=proc.time()-ptm
    stability$run.time=run_time
    solution_path=stability_run$x
    beta.full=apply(solution_path,1,max)
    beta.cv=beta.full
    beta.cv[beta.cv<0.6]=0
    
    tp=sum(beta.cv!=0 & beta0!=0)
    fp=sum(beta.cv!=0 & beta0==0)
    tn=sum(beta.cv==0 & beta0==0)
    fn=sum(beta.cv==0 & beta0!=0)
    if(tp+fp+tn+fn!=p) {warning("wrong tp/fp/tn/fn!")}
    if((tp+fp)==0) {mcc=0;ppv=0}
    if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
    tpr=tp/(tp+fn)
    stability$mcc=mcc
    stability$ppv=ppv
    stability$tpr=tpr
    
    pred=prediction(abs(beta.full),ifelse(beta0!=0,1,0))
    auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
    pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
    stability$pauc=pauc
    stability_list[[iter]]=stability
    
    
    ## for conference ##
    print("ciclos internos de 64:")
    print(iter)
    
  }
  
  ## Now, for the means and parameters ##
  
  r <- n/(s0*log(p-s0))
  
  G.scenario <- append(G.scenario, rep(cen, 7))
  G.method <- append(G.method, method_name)
  G.r <- append(G.r, rep(r, 7))
  G.n <- append(G.n, rep(n, 7))
  G.p <- append(G.p, rep(p, 7))
  G.s0 <- append(G.s0, rep(s0, 7))
  G.SNR <- append(G.SNR, rep(SNR, 7))
  G.corDesign <- append(G.corDesign, rep(corDesign, 7))
  G.rho <- append(G.rho, rep(rho, 7))
  G.pB <- append(G.pB, rep(pB, 7))
  G.s0B <- append(G.s0B, rep(s0B, 7))
  
  for(name in method_name) {assign(sprintf("%s_list_pauc",name),list(NA))}
  for(name in method_name) {assign(sprintf("%s_list_rmse",name),list(NA))}
  for(name in method_name) {assign(sprintf("%s_list_ppv",name),list(NA))}
  for(name in method_name) {assign(sprintf("%s_list_tpr",name),list(NA))}
  for(name in method_name) {assign(sprintf("%s_list_mcc",name),list(NA))}
  
  for(it in 1:3)
    {
    lasso_list_pauc[[1]] <- append(lasso_list_pauc[[1]], lasso_list[[it]]$pauc)
    lenet_list_pauc[[1]] <- append(lenet_list_pauc[[1]], lenet_list[[it]]$pauc)
    henet_list_pauc[[1]] <- append(henet_list_pauc[[1]], henet_list[[it]]$pauc)
    ridge_list_pauc[[1]] <- append(ridge_list_pauc[[1]], ridge_list[[it]]$pauc)
    dantzig_list_pauc[[1]] <- append(dantzig_list_pauc[[1]], NA)
    scad_list_pauc[[1]] <- append(scad_list_pauc[[1]], scad_list[[it]]$pauc)
    stability_list_pauc[[1]] <- append(stability_list_pauc[[1]], stability_list[[it]]$pauc)
  
    lasso_list_rmse[[1]] <- append(lasso_list_rmse[[1]], lasso_list[[it]]$rmse)
    lenet_list_rmse[[1]] <- append(lenet_list_rmse[[1]], lenet_list[[it]]$rmse)
    henet_list_rmse[[1]] <- append(henet_list_rmse[[1]], henet_list[[it]]$rmse)
    ridge_list_rmse[[1]] <- append(ridge_list_rmse[[1]], ridge_list[[it]]$rmse)
    dantzig_list_rmse[[1]] <- append(dantzig_list_rmse[[1]], NA)
    scad_list_rmse[[1]] <- append(scad_list_rmse[[1]], scad_list[[it]]$rmse)
    stability_list_rmse[[1]] <- append(stability_list_rmse[[1]], stability_list[[it]]$rmse)
    
    lasso_list_ppv[[1]] <- append(lasso_list_ppv[[1]],  lasso_list[[it]]$ppv)
    lenet_list_ppv[[1]] <- append(lenet_list_ppv[[1]],  lenet_list[[it]]$ppv)
    henet_list_ppv[[1]] <- append(henet_list_ppv[[1]],  henet_list[[it]]$ppv)
    ridge_list_ppv[[1]] <- append(ridge_list_ppv[[1]], ridge_list[[it]]$ppv)
    dantzig_list_ppv[[1]] <- append(dantzig_list_ppv[[1]], NA)
    scad_list_ppv[[1]] <- append(scad_list_ppv[[1]],  scad_list[[it]]$ppv)
    stability_list_ppv[[1]] <- append(stability_list_ppv[[1]],  stability_list[[it]]$ppv)

    lasso_list_mcc[[1]] <- append(lasso_list_mcc[[1]], lasso_list[[it]]$mcc)
    lenet_list_mcc[[1]] <- append(lenet_list_mcc[[1]], lenet_list[[it]]$mcc)
    henet_list_mcc[[1]] <- append(henet_list_mcc[[1]], henet_list[[it]]$mcc)
    ridge_list_mcc[[1]] <- append(ridge_list_mcc[[1]], ridge_list[[it]]$mcc)
    dantzig_list_mcc[[1]] <- append(dantzig_list_mcc[[1]], NA)[[1]]
    scad_list_mcc[[1]] <- append(scad_list_mcc[[1]], scad_list[[it]]$mcc)
    stability_list_mcc[[1]] <- append(stability_list_mcc[[1]], stability_list[[it]]$mcc)
    
    lasso_list_tpr[[1]] <- append(lasso_list_tpr[[1]], lasso_list[[it]]$tpr)
    lenet_list_tpr[[1]] <- append(lenet_list_tpr[[1]], lenet_list[[it]]$tpr)
    henet_list_tpr[[1]] <- append(henet_list_tpr[[1]], henet_list[[it]]$tpr)
    ridge_list_tpr[[1]] <- append(ridge_list_tpr[[1]], ridge_list[[it]]$tpr)
    dantzig_list_tpr[[1]] <- append(dantzig_list_tpr, NA)
    scad_list_tpr[[1]] <- append(scad_list_tpr[[1]], scad_list[[it]]$tpr)
    stability_list_tpr[[1]] <- append(stability_list_tpr[[1]], stability_list[[it]]$tpr)
  }
  
  m_lasso_list_rmse <- mean(lasso_list_rmse[[1]][2:4]); m_lenet_list_rmse <- mean(lenet_list_rmse[[1]][2:4])
  m_henet_list_rmse <- mean(henet_list_rmse[[1]][2:4]); m_ridge_list_rmse <- mean(ridge_list_rmse[[1]][2:4])
  m_dantzig_list_rmse <- NA; m_scad_list_rmse <- mean(scad_list_rmse[[1]][2:4])
  m_stability_list_rmse <- mean(stability_list_rmse[[1]][2:4])
  
  m_lasso_list_ppv <- mean(lasso_list_ppv[[1]][2:4]); m_lenet_list_ppv <- mean(lenet_list_ppv[[1]][2:4])
  m_henet_list_ppv <- mean(henet_list_ppv[[1]][2:4]); m_ridge_list_ppv <- mean(ridge_list_ppv[[1]][2:4])
  m_dantzig_list_ppv <- NA; m_scad_list_ppv <- mean(scad_list_ppv[[1]][2:4])
  m_stability_list_ppv <- mean(stability_list_ppv[[1]][2:4])
  
  m_lasso_list_mcc <- mean(lasso_list_mcc[[1]][2:4]); m_lenet_list_mcc <- mean(lenet_list_mcc[[1]][2:4])
  m_henet_list_mcc <- mean(henet_list_mcc[[1]][2:4]); m_ridge_list_mcc <- mean(ridge_list_mcc[[1]][2:4])
  m_dantzig_list_mcc <- NA; m_scad_list_mcc <- mean(scad_list_mcc[[1]][2:4])
  m_stability_list_mcc <- mean(stability_list_mcc[[1]][2:4])
  
  m_lasso_list_tpr <- mean(lasso_list_tpr[[1]][2:4]); m_lenet_list_tpr <- mean(lenet_list_tpr[[1]][2:4])
  m_henet_list_tpr <- mean(henet_list_tpr[[1]][2:4]); m_ridge_list_tpr <- mean(ridge_list_tpr[[1]][2:4])
  m_dantzig_list_tpr <- NA; m_scad_list_tpr <- mean(scad_list_tpr[[1]][2:4])
  m_stability_list_tpr <- mean(stability_list_tpr[[1]][2:4])
  
  m_lasso_list_pauc <- mean(lasso_list_pauc[[1]][2:4]); m_lenet_list_pauc <- mean(lenet_list_pauc[[1]][2:4])
  m_henet_list_pauc <- mean(henet_list_pauc[[1]][2:4]); m_ridge_list_pauc <- mean(ridge_list_pauc[[1]][2:4])
  m_dantzig_list_pauc <- NA; m_scad_list_pauc <- mean(scad_list_pauc[[1]][2:4])
  m_stability_list_pauc <- mean(stability_list_pauc[[1]][2:4])
  
  s_lasso_list_rmse <- sd(lasso_list_rmse[[1]][2:4]); s_lenet_list_rmse <- sd(lenet_list_rmse[[1]][2:4])
  s_henet_list_rmse <- sd(henet_list_rmse[[1]][2:4]); s_ridge_list_rmse <- sd(ridge_list_rmse[[1]][2:4])
  s_dantzig_list_rmse <- NA; s_scad_list_rmse <- sd(scad_list_rmse[[1]][2:4])
  s_stability_list_rmse <- sd(stability_list_rmse[[1]][2:4])
  
  s_lasso_list_ppv <- sd(lasso_list_ppv[[1]][2:4]); s_lenet_list_ppv <- sd(lenet_list_ppv[[1]][2:4])
  s_henet_list_ppv <- sd(henet_list_ppv[[1]][2:4]); s_ridge_list_ppv <- sd(ridge_list_ppv[[1]][2:4])
  s_dantzig_list_ppv <- NA; s_scad_list_ppv <- sd(scad_list_ppv[[1]][2:4])
  s_stability_list_ppv <- sd(stability_list_ppv[[1]][2:4])
  
  s_lasso_list_mcc <- sd(lasso_list_mcc[[1]][2:4]); s_lenet_list_mcc <- sd(lenet_list_mcc[[1]][2:4])
  s_henet_list_mcc <- sd(henet_list_mcc[[1]][2:4]); s_ridge_list_mcc <- sd(ridge_list_mcc[[1]][2:4])
  s_dantzig_list_mcc <- NA; s_scad_list_mcc <- sd(scad_list_mcc[[1]][2:4])
  s_stability_list_mcc <- sd(stability_list_mcc[[1]][2:4])
  
  s_lasso_list_tpr <- sd(lasso_list_tpr[[1]][2:4]); s_lenet_list_tpr <- sd(lenet_list_tpr[[1]][2:4])
  s_henet_list_tpr <- sd(henet_list_tpr[[1]][2:4]); s_ridge_list_tpr <- sd(ridge_list_tpr[[1]][2:4])
  s_dantzig_list_tpr <- NA; s_scad_list_tpr <- sd(scad_list_tpr[[1]][2:4])
  s_stability_list_tpr <- sd(stability_list_tpr[[1]][2:4])
  
  s_lasso_list_pauc <- sd(lasso_list_pauc[[1]][2:4]); s_lenet_list_pauc <- sd(lenet_list_pauc[[1]][2:4])
  s_henet_list_pauc <- sd(henet_list_pauc[[1]][2:4]); s_ridge_list_pauc <- sd(ridge_list_pauc[[1]][2:4])
  s_dantzig_list_pauc <- NA; s_scad_list_pauc <- sd(scad_list_pauc[[1]][2:4])
  s_stability_list_pauc <- sd(stability_list_pauc[[1]][2:4])
  
  G.mean_pauc <- append(G.mean_pauc, c(m_lasso_list_pauc, m_lenet_list_pauc, m_henet_list_pauc,
                                       m_ridge_list_pauc, m_dantzig_list_pauc, m_scad_list_pauc,
                                       m_stability_list_pauc))
  G.sd_pauc <- append(G.sd_pauc,  c(s_lasso_list_pauc, s_lenet_list_pauc, s_henet_list_pauc,
                                    s_ridge_list_pauc, s_dantzig_list_pauc, s_scad_list_pauc,
                                    s_stability_list_pauc))
  G.mean_rmse <- append(G.mean_rmse,  c(m_lasso_list_rmse, m_lenet_list_rmse, m_henet_list_rmse,
                                        m_ridge_list_rmse, m_dantzig_list_rmse, m_scad_list_rmse,
                                        m_stability_list_rmse))
  G.sd_rmse <- append(G.sd_rmse,  c(s_lasso_list_rmse, s_lenet_list_rmse, s_henet_list_rmse,
                                    s_ridge_list_rmse, s_dantzig_list_rmse, s_scad_list_rmse,
                                    s_stability_list_rmse))
  G.mean_tpr <- append(G.mean_tpr,  c(m_lasso_list_tpr, m_lenet_list_tpr, m_henet_list_tpr,
                                      m_ridge_list_tpr, m_dantzig_list_tpr, m_scad_list_tpr,
                                      m_stability_list_tpr))
  G.sd_tpr <- append(G.sd_tpr,  c(s_lasso_list_tpr, s_lenet_list_tpr, s_henet_list_tpr,
                                  s_ridge_list_tpr, s_dantzig_list_tpr, s_scad_list_tpr,
                                  s_stability_list_tpr))
  G.mean_ppv <- append(G.mean_ppv,  c(m_lasso_list_ppv, m_lenet_list_ppv, m_henet_list_ppv,
                                      m_ridge_list_ppv, m_dantzig_list_ppv, m_scad_list_ppv,
                                      m_stability_list_ppv))
  G.sd_ppv <- append(G.sd_ppv,  c(s_lasso_list_ppv, s_lenet_list_ppv, s_henet_list_ppv,
                                  s_ridge_list_ppv, s_dantzig_list_ppv, s_scad_list_ppv,
                                  s_stability_list_ppv))
  G.mean_mcc <- append(G.mean_mcc,  c(m_lasso_list_mcc, m_lenet_list_mcc, m_henet_list_mcc,
                                      m_ridge_list_mcc, m_dantzig_list_mcc, m_scad_list_mcc,
                                      m_stability_list_mcc))
  G.sd_mcc <- append(G.sd_mcc,  c(s_lasso_list_mcc, s_lenet_list_mcc, s_henet_list_mcc,
                                  s_ridge_list_mcc, s_dantzig_list_mcc, s_scad_list_mcc,
                                  s_stability_list_mcc))
  
  
  if(cen == 1){
    G.mean_pauc <- G.mean_pauc[-1]; G.sd_pauc <- G.sd_pauc[-1];G.mean_rmse <- G.mean_rmse[-1]
    G.sd_rmse <- G.sd_rmse[-1]; G.mean_tpr <- G.mean_tpr[-1]; G.sd_tpr <- G.sd_tpr[-1]
    G.mean_ppv <- G.mean_ppv[-1]; G.sd_ppv <- G.sd_ppv[-1]; G.mean_mcc <- G.mean_mcc[-1]
    G.sd_mcc <- G.sd_mcc[-1] 
  }
  
  print("Ciclo externo:")
  print(cen)
  
}

pauc_final <- data.frame(G.scenario, G.method, G.r, G.n, G.p, G.s0, G.SNR, G.corDesign,
                         G.rho, G.pB, G.s0B, G.mean_pauc, G.sd_pauc)
rmse_final <- data.frame(G.scenario, G.method, G.r, G.n, G.p, G.s0, G.SNR, G.corDesign,
                        G.rho, G.pB, G.s0B, G.mean_rmse, G.sd_rmse)
tpr_final <- data.frame(G.scenario, G.method, G.r, G.n, G.p, G.s0, G.SNR, G.corDesign,
                       G.rho, G.pB, G.s0B, G.mean_tpr, G.sd_tpr)
ppv_final <- data.frame(G.scenario, G.method, G.r, G.n, G.p, G.s0, G.SNR, G.corDesign,
                       G.rho, G.pB, G.s0B, G.mean_ppv, G.sd_ppv)
mcc_final <- data.frame(G.scenario, G.method, G.r, G.n, G.p, G.s0, G.SNR, G.corDesign,
                        G.rho, G.pB, G.s0B, G.mean_mcc, G.sd_mcc)

scores_local_syn <- list()
scores_local_syn <- append(scores_local_syn, list(pauc_final, rmse_final, tpr_final, ppv_final, mcc_final))
names(scores_local_syn) <- c('pauc', 'rmse', 'tpr', 'ppv', 'mcc')

