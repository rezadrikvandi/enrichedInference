# Enriched estimation and inference for high dimensional data

#---------------------------------------------------------------------------------------------------
#' This function calculates enriched estimates and inference for a high dimensional regression model
#' @param X                  an nxp matrix of covariates (i.e., n observations and p covariates)
#' @param Y                  a vector of size n containing response values
#' @param k                  number of principal components for the supervised SVD (default is k=1)
#' @param intercept          intercept for the regression model (default is TRUE)
#' @param standardize        standardise cavariates (default is FALSE)
#' @param CVlasso            whether to use CV for lasso tuning parameter selection (default is FALSE)
#' @param Xcorr              use of either "marginal" or "partial" correlation in the supervised SVD
#' @param sig.level          significance level for enriched inferences (default is 0.05)
#' @param Bonf.correction    Bonferroni correction for enriched multiple testing (default is FALSE)
#' @param threshold          threshold for the minimum correlation in the supervised SVD (default is 0.1)
#'
#' @return    enriched estimates of regression parameters for selected covariates and their standard errors,
#'            enriched confidence intervals, p-values for enriched tests,
#'            enriched estimate of error variance, estimate of thetak, k, and negative 2log-likelihood of the enriched selected model
#' @export

enrichedInference <- function(X, Y, k=1, intercept=TRUE, standardize=FALSE, CVlasso=FALSE, Xcorr=c("marginal","partial"), threshold=0.1, sig.level=0.05, Bonf.correction=FALSE)
{
  Y.train  <- as.vector(Y);
  X.train <- as.matrix(X);
  nn <- nrow(X.train);
  p <- ncol(X.train);
  rownames(X.train) <- NULL;
  colnames(X.train) <- NULL;
  dimnames(X.train) <- list(rownames(X.train, do.NULL=TRUE),colnames(X.train, do.NULL=FALSE, prefix = "X"));
  if(standardize==TRUE)
  {
    X.train <- scale(X.train,center=TRUE,scale=TRUE);
    M0 <- glmnet(X.train,Y.train, family="gaussian", alpha=1, standardize=FALSE);
    CVL <- cv.glmnet(X.train,Y.train, family="gaussian", alpha=1, standardize=FALSE);
    LambdaMinCV <- CVL$lambda.min;
    betahatCV <- coef(M0, x=X.train,y=Y.train, s=LambdaMinCV, standardize=FALSE, exact=TRUE);
  } else if(standardize==FALSE){
    M0 <- glmnet(X.train,Y.train, family="gaussian", alpha=1);
    CVL <- cv.glmnet(X.train,Y.train, family="gaussian", alpha=1);
    LambdaMinCV <- CVL$lambda.min;
    betahatCV <- coef(M0, x=X.train,y=Y.train, s=LambdaMinCV, exact=TRUE);
  }
  bCVnew <- apply(betahatCV, 1, function(row) all(row !=0));
  estimateCV <- betahatCV[bCVnew,];
  estimateCV <- cbind(estimateCV);
  Index_XsCV <- row.names(estimateCV);
  Index_XsCV <- c(Index_XsCV);
  Index_XsCV <- gsub("[a-zA-Z ]", "", Index_XsCV);
  Index_XsCV <- Index_XsCV[-1];
  Index_XsCV <- as.numeric(Index_XsCV);
  q <- length(Index_XsCV);
  Xs <- X.train[,Index_XsCV];
  dimnames(Xs) <- list(rownames(Xs, do.NULL=TRUE),colnames(Xs, do.NULL=FALSE, prefix = "Xs"));
  Int <- rep(1,nn);
  Xs <- cbind(Int,Xs);
  sigma2CV <- t(Y.train-Xs%*%as.vector(estimateCV))%*%((Y.train-Xs%*%as.vector(estimateCV)))/(nn-(q+1));
  sigma2CV <- max(0,sigma2CV);
  sigma2CV[sigma2CV==Inf] <- 0;
  I_n <- diag(nn,x=rep(1,nn));
  CovEps <- c(sigma2CV)*I_n;
  if(semidefiniteness(CovEps))
  {
    CovEps <- CovEps;
  } else{
    eigCovEps <- eigen(CovEps, symmetric=TRUE);
    eigCovEps$values <- pmax(0, eigCovEps$values);
    CovEps <- eigCovEps$vectors%*%diag(eigCovEps$values)%*%t(eigCovEps$vectors);
  }
  eps <- rmvnorm(10000,mean=rep(0,nn),sigma=CovEps);
  if(CVlasso==FALSE & standardize==TRUE)
  {
    Xeps <- t(X.train)%*%t(eps);
    infitynorm_Xeps <- apply(Xeps,2,max);
    Lambda_Lasso <- 2*mean(infitynorm_Xeps)/nn;
    M1 <- glmnet(X.train,Y.train, alpha=1, family="gaussian", standardize=FALSE);
    betahat1 <- coef(M1, x=X.train,y=Y.train, s=Lambda_Lasso, standardize=FALSE, exact=TRUE);
  } else if(CVlasso==FALSE & standardize==FALSE){
    X.train2 <- scale(X.train,center=TRUE,scale=FALSE);
    Xeps <- t(X.train2)%*%t(eps);
    infitynorm_Xeps <- apply(Xeps,2,max);
    Lambda_Lasso <- 2*mean(infitynorm_Xeps)/nn;
    M1 <- glmnet(X.train,Y.train, alpha=1, family="gaussian");
    betahat1 <- coef(M1, x=X.train,y=Y.train, s=Lambda_Lasso, exact=TRUE);
  } else if(CVlasso==TRUE){
    Lambda_Lasso <- LambdaMinCV;
    betahat1 <- betahatCV;
  }
  betahat1 <- as.matrix(betahat1);
  b1new <- apply(betahat1, 1, function(row) all(row!=0));
  estimate1 <- betahat1[b1new,];
  leng <- length(estimate1);
  if(leng<2)
  {
    stop("For this dataset, the lasso did not select any covariates, so the enriched method is not applied.", call.=FALSE)
  }
  Lasso_estimates <- as.matrix(estimate1);
  Index_Xs <- row.names(Lasso_estimates);
  Index_Xs <- c(Index_Xs);
  Index_Xs <- gsub("[a-zA-Z ]", "", Index_Xs);
  Index_Xs <- Index_Xs[-1];
  Index_Xs <- as.numeric(Index_Xs);
  Signs <- as.vector(ifelse(Lasso_estimates>0,1,-1));
  Index_of_Xs <- Index_Xs;
  q <- length(Index_of_Xs);
  length_Index_Xs <- q;
  Xs <- X.train[,Index_of_Xs];
  Xu <- X.train[,-Index_of_Xs];

  if(ncol(Xu)==0 | is.null(ncol(Xu))){
    print(Lasso_estimates)
    stop("For this dataset, all the covariates are selected by lasso and there is no unselected covariate. So, the lasso estimates are returned here.", call.=FALSE)
  }

  Xs <- matrix(c(Xs),nn,length_Index_Xs);
  dimnames(Xs) <- list(rownames(Xs, do.NULL=TRUE),colnames(Xs, do.NULL=FALSE, prefix = "Xs"));

  remove_Xs <- p;
  if(q>nn)
  {
    CXs <- apply(Xs, 2, function(a) (t(a)%*%Y.train/sqrt(t(a)%*%a)))
    rank_CXs <- rank(-abs(as.vector(CXs)))
    rank_est <- rank(-abs(as.vector(Lasso_estimates[-1])))
    cp <- ifelse(p<5000,0.7,0.3)
    remove_Xs <- which(rank_CXs>q*cp & rank_est>q*cp)
    if(!identical(remove_Xs,integer(0)))
    {
      Xu <- cbind(Xs[,remove_Xs],Xu)
      Xs <- Xs[,-remove_Xs]
      q <- ncol(Xs)
      Signs <- c(Signs[1],Signs[-1][-remove_Xs])
      Index_Xs <- Index_of_Xs[-remove_Xs]
    } else if(identical(remove_Xs,integer(0))){remove_Xs <- ncol(X.train)}
  }

  if(is.null(q)){
    q <- 1
    Xs <- matrix(c(Xs),nn,1)
    dimnames(Xs) <- list(rownames(Xs, do.NULL = TRUE),colnames(Xs, do.NULL = FALSE, prefix = "Xs"));
  }

  if(intercept==TRUE)
  {
    Int <- rep(1,nn);
    Xs <- cbind(Int,Xs);
    p <- ncol(X.train)+1;
  } else{
    Signs <- Signs[-1];
    p <- ncol(X.train);
  }

  if(!is.null(colnames(X)))
  {
    selectedX <- colnames(X)[Index_Xs];
  } else{
    selectedX <- Index_Xs;
  }

  if(Xcorr=="partial")
  {
    if(ncol(Xu)>5000)
    {
      Xu_sd <- scale(Xu, center = TRUE,scale = TRUE)
      Y.train_sd <- as.matrix(Y.train,ncol=1)
      Y.train_sd <- scale(Y.train_sd, center = TRUE,scale = TRUE)
      CXu <- apply(Xu_sd, 2, function(a) (t(a)%*%Y.train_sd/(sqrt(t(a)%*%a)*sqrt(t(Y.train_sd)%*%Y.train_sd))));
      CXu <- abs(CXu)
      remove_Xu <- which(CXu<threshold)
      Xu <- Xu[,-remove_Xu]
    }
    pcorr <- suppressWarnings(try(pcor(cbind(Y.train,Xu)),silent=TRUE));
    Xu <- Xu[,-which(c(pcorr$estimate[1,])[-1]<threshold)];
  }
  if(Xcorr=="marginal")
  {
    Xu_sd <- scale(Xu, center = TRUE,scale = TRUE)
    Y.train_sd <- as.matrix(Y.train,ncol=1)
    Y.train_sd <- scale(Y.train_sd, center = TRUE,scale = TRUE)
    CXu <- apply(Xu_sd, 2, function(a) (t(a)%*%Y.train_sd/(sqrt(t(a)%*%a)*sqrt(t(Y.train_sd)%*%Y.train_sd))));
    CXu <- abs(CXu)
    remove_Xu <- which(CXu<threshold)
    Xu <- Xu[,-remove_Xu]
  }
  q <- ncol(Xs);
  pp <- ncol(Xu)+q;

  Xu_decomposition=svd(Xu,nv=ncol(Xu));

  if(ncol(Xu)-length(Xu_decomposition$d)!=0)
  {
    add1=matrix(0,length(Xu_decomposition$d),ncol(Xu)-length(Xu_decomposition$d));
    add2=matrix(0,ncol(Xu)-length(Xu_decomposition$d),ncol(Xu));
    U=cbind(Xu_decomposition$u,add1);
    D=rbind(cbind(diag(Xu_decomposition$d),add1),add2);
    P=U%*%D;
    V=Xu_decomposition$v;
  } else{
    U=Xu_decomposition$u;
    D=diag(Xu_decomposition$d);
    P=U%*%D;
    V=Xu_decomposition$v;
  }

  if(is.null(k))
  {
    k <- 1;
  } else{
    k <- floor(k);
  }

  if(k<=0)
  {
    stop("k cannot be set to 0 or negative. Use k=NULL, the default value of k, or a positive integer k", call.=FALSE)
  }

  if(k>=1 & k<=ncol(U))
  {
    U_k <- U[,1:k];
    D_k <- D[1:k,1:k];
    V_k <- V[1:k,1:k];
    P_k <- P[,1:k];
  } else{
    U_k <- U;
    D_k <- D;
    V_k <- V;
    P_k <- P;
  }

  if(k==1){P_k <- as.matrix(P_k)}
  dimnames(P_k) <- list(rownames(P_k, do.NULL = TRUE),colnames(P_k, do.NULL = FALSE, prefix = "PC"));

  profile_likelihood <- function(par)
  {
    lambdaa <- 0.5/(exp(par)+1);
    I_n=diag(nn,x=rep(1,nn));
    Sigma_lambdaa <- (1-lambdaa)*I_n + lambdaa*U_k%*%t(U_k);
    Z <- Xs;
    H <- Z%*%ginv(t(Z)%*%ginv(Sigma_lambdaa)%*%Z)%*%t(Z)%*%ginv(Sigma_lambdaa);
    e <- (I_n-H)%*%Y.train;
    l_p=(-1/2)*(log(det(Sigma_lambdaa)))-(nn/2)*(log(t(e)%*%ginv(Sigma_lambdaa)%*%e));
    return(c(-l_p));
  }

  par_hat_unconditionalEPoS <- suppressWarnings(try(nlminb(log(0.5/0.1-1),profile_likelihood,control=list(iter.max=1000))$par,silent=TRUE));
  if(typeof(par_hat_unconditionalEPoS)=="character")
  {
    Z <- Xs;
    betahat_s_marginal_unconditionalEPoS <- c(ginv(t(Z)%*%Z)%*%t(Z)%*%Y.train);
    lambda_initial <- 0.1;
    sigma2_initial <- 1;
    beta_s_initial <- betahat_s_marginal_unconditionalEPoS;
  } else{
    if(is.na(par_hat_unconditionalEPoS) | exp(par_hat_unconditionalEPoS)==-Inf)
    {
      Z <- Xs;
      betahat_s_marginal_unconditionalEPoS <- c(ginv(t(Z)%*%Z)%*%t(Z)%*%Y.train);
      lambda_initial <- 0.1;
      sigma2_initial <- 1;
      beta_s_initial <- betahat_s_marginal_unconditionalEPoS;
    } else{
      lambda_hat_unconditionalEPoS <- 0.5/(exp(par_hat_unconditionalEPoS)+1);
      I_n=diag(nn,x=rep(1,nn));
      Sigma_lambda_hat_unconditionalEPoS <- (1-lambda_hat_unconditionalEPoS)*I_n + lambda_hat_unconditionalEPoS*U_k%*%t(U_k);
      Z <- Xs;
      betahat_s_marginal_unconditionalEPoS <- ginv(t(Z)%*%ginv(Sigma_lambda_hat_unconditionalEPoS)%*%Z)%*%t(Z)%*%ginv(Sigma_lambda_hat_unconditionalEPoS)%*%Y.train;
      sigma2_hat_marginal_unconditionalEPoS <- c((1/nn)*t(Y.train-Xs%*%betahat_s_marginal_unconditionalEPoS)%*%(Y.train-Xs%*%betahat_s_marginal_unconditionalEPoS));
      beta_s_initial  <- c(betahat_s_marginal_unconditionalEPoS);
      sigma2_initial  <- c(sigma2_hat_marginal_unconditionalEPoS);
      lambda_initial <- c(lambda_hat_unconditionalEPoS);
    }}

  I_n <- diag(nn,x=rep(1,nn));
  Ip_q <- matrix(rep(1,pp-q),pp-q,1);
  A1 <- t(Xu)%*%(I_n-Xs%*%ginv(t(Xs)%*%Xs)%*%t(Xs));
  A1 <- (1/(nn*Lambda_Lasso))*A1;
  A2 <- -A1;
  A3 <- -diag(Signs,nrow=length(Signs),ncol=length(Signs))%*%ginv(t(Xs)%*%Xs)%*%t(Xs);
  A <- rbind(A1,A2,A3);
  B1 <- Ip_q - t(Xu)%*%ginv(Xs%*%t(Xs))%*%Xs%*%Signs;
  B2 <- Ip_q + t(Xu)%*%ginv(Xs%*%t(Xs))%*%Xs%*%Signs;
  B3 <- -diag(Signs,nrow=length(Signs),ncol=length(Signs))%*%ginv(t(Xs)%*%Xs)%*%Signs;
  B3 <- (nn*Lambda_Lasso)*B3;
  B <- rbind(B1,B2,B3);
  KKT1 <- max(A1%*%Y.train-B1);
  KKT2 <- max(A2%*%Y.train-B2);
  KKT3 <- max(A3%*%Y.train-B3);
  A <- A3;
  B <- B3;
  max(A%*%Y.train-B);

  betaas <- beta_s_initial;
  sigmaa2 <- sigma2_initial;
  lambdaa <- lambda_initial;
  I_n=diag(nn,x=rep(1,nn));
  Sigma_lambdaa <- (1-lambdaa)*I_n + lambdaa*U_k%*%t(U_k);
  mean_y <- Xs%*%betaas;
  cov_y <- (sigmaa2)*Sigma_lambdaa;
  n_r <- 100000;
  y_r <- rmvnorm(n_r,mean=c(mean_y),sigma=cov_y);
  TF <- apply(y_r, 1, function(t) (max(A%*%t-B)<=0));
  E <- length(TF[TF==TRUE]);
  y_final <- matrix(y_r[TF,],E,nn);
  y_final2 <- t(as.matrix(y_final));
  if(ncol(y_final2)>1)
  {
    lo_lim <- apply(y_final2, 1, FUN=min);
    up_lim <- apply(y_final2, 1, FUN=max);
    ylim <-cbind(lo_lim,up_lim);
  }else{ylim <-cbind(rep(-Inf,nn),rep(Inf,nn));}

  conditional_EPoS_likelihood <- function(par)
  {
    lambdaa <- 0.5/(exp(par[1])+1); sigmaa2 <- exp(par[2]); betas <- par[3:(q+2)];
    I_n=diag(nn,x=rep(1,nn));
    Sigma_lambdaa <- (1-lambdaa)*I_n + lambdaa*U_k%*%t(U_k);
    mean_y <- Xs%*%betas;
    cov_y <- (sigmaa2)*Sigma_lambdaa;
    l_p <- dmvnorm(x=Y.train,mean=c(mean_y),sigma=cov_y,log=TRUE)-log(pmvnorm(lower=c(ylim[,1]),upper=c(ylim[,2]),mean=c(mean_y),sigma=cov_y,maxpts=10000)[1]);
    return(c(-l_p));
  }

  par_hat <- try(nlminb(c(log(0.5/lambda_initial-1),log(sigma2_initial),beta_s_initial),conditional_EPoS_likelihood,control=list(iter.max=100))$par,silent=TRUE);
  if(typeof(par_hat)=="character")
  {
    npar <- length(c(beta_s_initial))+2;
    par_hat <- suppressWarnings(try(summary(optimx(par=c(log(1/lambda_initial-1),log(sigma2_initial),beta_s_initial),fn=conditional_EPoS_likelihood,control=list(all.methods=TRUE,kkt=FALSE,maxit=100)), order=value)[1,1:npar],silent=TRUE));
  }

  if(typeof(par_hat)=="character")
  {
    cat("Use another value of k instead of the default k, or find an optimal k using the Cross Validation function in the R package.");
  } else{
    par_hat <- as.vector(unlist(par_hat));
    lambda_hat <- 0.5/(exp(par_hat[1])+1);
    sigma2_hat_marginal <- exp(par_hat[2]);
    betahat_s_marginal <- par_hat[-c(1:2)];

    neg2loglik <- 2*(conditional_EPoS_likelihood(par_hat));

    conditional_EPoS_likelihood2 <- function(par)
    {
      lambdaa <- 0.5/(exp(par_hat[1])+1); sigmaa2 <- exp(par_hat[2]); betas <- par[c(1:q)];
      I_n=diag(nn,x=rep(1,nn));
      Sigma_lambdaa <- (1-lambdaa)*I_n + lambdaa*U_k%*%t(U_k);
      mean_y <- Xs%*%betas;
      cov_y <- (sigmaa2)*Sigma_lambdaa;
      l_p <- dmvnorm(x=Y.train,mean=c(mean_y),sigma=cov_y,log=TRUE)-log(pmvnorm(lower=c(ylim[,1]),upper=c(ylim[,2]),mean=c(mean_y),sigma=cov_y,maxpts=10000)[1]);
      return(c(-l_p));
    }

    hessian_negative <- hessian(x=par_hat[-c(1,2)], func=conditional_EPoS_likelihood2);
    if(!is.na(sum(hessian_negative)) & sum(hessian_negative)!=Inf & sum(hessian_negative)!=-Inf)
    {
      var_betahat_s_marginal <- ginv(hessian_negative);
    } else
    {
      Sigma_lambda_hat <- (1-lambda_hat)*I_n + lambda_hat*U_k%*%t(U_k);
      Z <- Xs;
      var_betahat_s_marginal <- sigma2_hat_marginal*ginv(t(Z)%*%ginv(Sigma_lambda_hat)%*%Z);
    }
    var_betahat_s_marginal <- round(var_betahat_s_marginal, digits=6);
    if(semidefiniteness(var_betahat_s_marginal))
    {
      var_betahat_s_marginal <- var_betahat_s_marginal;
    } else{
      eig <- suppressWarnings(try(eigen(var_betahat_s_marginal, symmetric=TRUE),silent=TRUE));
      eig$values <- pmax(0, eig$values);
      var_betahat_s_marginal <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors);
    }

    if(Bonf.correction==FALSE)
    {
      alpha <- sig.level;
    } else{
      alpha <- ifelse(intercept==TRUE,sig.level/length(betahat_s_marginal[-1]),sig.level/length(betahat_s_marginal));
    }

    CI <- cbind(betahat_s_marginal-qnorm(1-alpha/2,0,1)*sqrt(diag(var_betahat_s_marginal)),betahat_s_marginal+qnorm(1-alpha/2,0,1)*sqrt(diag(var_betahat_s_marginal)));

    TestStat <- numeric(length(betahat_s_marginal));
    pvalue <- numeric(length(betahat_s_marginal));
    SE <- numeric(length(betahat_s_marginal));
    for(j in 1:length(betahat_s_marginal))
    {
      SE[j] <- sqrt(diag(var_betahat_s_marginal)[j]);
      TestStat[j] <- (betahat_s_marginal[j])/SE[j];
      pvalue[j] <- 2*(1-pnorm(abs(TestStat[j])));
    }
    if(Bonf.correction==TRUE & intercept==TRUE)
    {
      pvalue <- pvalue*length(betahat_s_marginal[-1]);
      pvalue[pvalue >= 1] <- 1;
    } else if(Bonf.correction==TRUE & intercept==FALSE){
      pvalue <- pvalue*length(betahat_s_marginal);
      pvalue[pvalue >= 1] <- 1;
    }
    betahat_s_marginal <- round(betahat_s_marginal, digits=4);
    SE <- round(SE, digits=4);
    CI <- round(CI, digits=4);
    TestStat <- round(TestStat, digits=4);
    pvalue <- signif(pvalue,5);
    stars <- add.significance.stars(pvalue, cutoffs=c(0.05, 0.01, 0.001));
    stars <- str_replace_all(stars, " ", "");
    pvalue <- interaction(pvalue,stars,sep = "",drop=TRUE);

    I_n <- diag(nn,x=rep(1,nn));
    Sigma_lambda_hat <- (1-lambda_hat)*I_n + lambda_hat*U_k%*%t(U_k);
    theta_hat <- lambda_hat*(ginv(D_k)%*%ginv(D_k)%*%t(P_k)%*%ginv(Sigma_lambda_hat))%*%(Y.train-Xs%*%betahat_s_marginal);
    sigma2_hat_marginal <- (1-lambda_hat)*sigma2_hat_marginal;

    if(intercept==TRUE)
    {
      results <- data.frame(c("Intercept",selectedX),c(betahat_s_marginal),c(SE),c(CI[,1]),c(CI[,2]),c(pvalue))
      colnames(results) <- c("Coefficients", "Estimate", "Std.Error","Lower 95% CI","Upper 95% CI","p-value");
      orig.names <- names(results);
      name.width <- max(sapply(names(results), nchar));
      names(results) <- format(names(results), width = name.width, justify = "right");
      output <- format(results, width = name.width, justify = "left");
      ColumnSpaced <- data.frame("p-value" = paste("  ", output$`     p-value`))
      output$`     p-value` <- ColumnSpaced$p.value;
      print(output,row.names=FALSE);
    } else{
      results <- data.frame(c(selectedX),c(betahat_s_marginal),c(SE),c(CI[,1]),c(CI[,2]),c(pvalue))
      colnames(results) <- c("Coefficients", "Estimate", "Std.Error","Lower 95% CI","Upper 95% CI","p-value");
      orig.names <- names(results);
      name.width <- max(sapply(names(results), nchar));
      names(results) <- format(names(results), width = name.width, justify = "right");
      output <- format(results, width = name.width, justify = "left");
      ColumnSpaced <- data.frame("p-value" = paste("  ", output$`     p-value`))
      output$`     p-value` <- ColumnSpaced$p.value;
      print(output,row.names=FALSE);
    }

    cat("\n");
    cat("The estimate of error variance, sigma2, is:", sigma2_hat_marginal);
    cat("\n");
    cat("The estimate of parameter(s) theta is:", theta_hat);
    cat("\n");
    cat("The value of k is:", k);
    cat("\n");
    cat("-2loglikelihood is:", neg2loglik);
    cat("\n");
    cat("\n");
  }

  if(intercept==TRUE)
  {
    return(list(Selected_Covariates=c("Intercept",selectedX),Estimate=betahat_s_marginal,Std.Error=SE,sigma2_hat=sigma2_hat_marginal,theta_hat=theta_hat,k=k,neg2loglik=neg2loglik));
  } else{
    return(list(Selected_Covariates=c(selectedX),Estimate=betahat_s_marginal,Std.Error=SE,sigma2_hat=sigma2_hat_marginal,theta_hat=theta_hat,k=k,neg2loglik=neg2loglik));
  }
}
