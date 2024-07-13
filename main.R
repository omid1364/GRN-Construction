 #library(pcalg)
## define parameters
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("Rgraphviz", version = "3.8")
library(glmnet)
	library("MASS") 
    set.seed(98989)
    amel<-read.table("C://MyFolders//D//resent paper//Adjusted_GRN//R code//S3_File//Gene expressions//amel.txt");
    ana<-read.table("C://MyFolders//D//resent paper//Adjusted_GRN//R code//S3_File//Gene expressions//ana.txt");
    per<-read.table("C://MyFolders//D//resent paper//Adjusted_GRN//R code//S3_File//Gene expressions//per.txt");
    pse<-read.table("C://MyFolders//D//resent paper//Adjusted_GRN//R code//S3_File//Gene expressions//pse.txt");
    sim<-read.table("C://MyFolders//D//resent paper//Adjusted_GRN//R code//S3_File//Gene expressions//sim.txt");
    vir<-read.table("C://MyFolders//D//resent paper//Adjusted_GRN//R code//S3_File//Gene expressions//vir.txt");
####################################   
    X<-t(amel)
    Y=t(per)
###### bootsrap sampling
   boot_y<-sample(1:dim(Y)[1],1000-dim(Y)[1], rep=T)
   boot_x<-sample(1:dim(X)[1],1000-dim(X)[1], rep=T)
   my_n1_y <- 1000-dim(Y)[1]                     # Specify sample size for y
   my_n1_x <- 1000-dim(X)[1]                     # Specify sample size for x
   my_mu1 <- rep(0,2049)               # Specify the means of the variables
   my_Sigma1 <- 0.001*diag(2049)        # Specify the covariance matrix of the variables
   noise_y<-mvrnorm(n = my_n1_y, mu = my_mu1, Sigma = my_Sigma1)
   noise_x<-mvrnorm(n = my_n1_x, mu = my_mu1, Sigma = my_Sigma1)               
 #####################################################
      
  YY<-Y[boot_y,]+noise_y 
    Y1<-rbind(Y,YY)
    #X_initial=t(read.csv('X.csv',header=F))
    #ss<-rowSums(X_initial)
    #X1<-X_initial[rowSums(X_initial)>44,]
   
    XX<-X[boot_x,]+noise_x 
    X1<-rbind(X,XX)


    n=dim(Y1)[1]
    p=dim(Y1)[2]
    q=dim(X1)[2]
                   
    X11=as.matrix(X1)
    Y11=as.matrix(Y1)
	
     X11=cbind(1,X11)

     set.seed(112)
     train_rows <- sample(1:n, .9*n)
     x.train <- X11[train_rows, ]
     x.test <- X11[-train_rows, ]

      mse_lasso <- c(rep(0,p))
      mse_ridge <- c(rep(0,p))
      mse_elanet<- c(rep(0,p))

           
	for(i in 1:p)
	{
       Yl=Y11[,i]
       y.train <- Yl[train_rows]
       y.test <- Yl[-train_rows]
      # Fit models:
      #lambda.grid = 10^seq(5, -2, length=100)  # lambda=lambda.grid
      fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1)
      fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0)
      fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5)
      # 10-fold Cross validation for each alpha = 0,0.5,1.0
      fit1 <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=1, 
                          family="gaussian")
      fit0 <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=0,
                          family="gaussian")
      fit5 <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=.5,
                          family="gaussian")

      #ridge.bestlam = fit1$lambda.min
      #ridge.lam1se = fit1$lambda.1se
      #lasso.bestlam = fit0$lambda.min
      #lasso.lam1se = fit0$lambda.1se
      #elastic.bestlam = fit5$lambda.min
      #elastic.lam1se = fit5$lambda.1se

      #ridge.mod.best <- glmnet(x.train, y.train, family="gaussian", alpha=0, lambda=ridge.bestlam)
      #ridge.mod.1se <- glmnet(x.train, y.train, family="gaussian", alpha=0, lambda=ridge.lam1se)
      #lasso.mod.best <- glmnet(x.train, y.train, family="gaussian", alpha=1, lambda=ridge.bestlam)
      #lasso.mod.1se <- glmnet(x.train, y.train, family="gaussian", alpha=1, lambda=ridge.lam1se)
      #elastic.mod.best <- glmnet(x.train, y.train, family="gaussian", alpha=.5, lambda=ridge.bestlam)
      #elastic.mod.1se <- glmnet(x.train, y.train, family="gaussian", alpha=.5, lambda=ridge.lam1se)


       #coef(ridge.mod.best)
       #coef(ridge.mod.1se)
       #coef(lasso.mod.best)
       #coef(lasso.mod.1se)
       #coef(elastic.mod.best)
       #coef(elastic.mod.1se)


      # Plot solution paths:
       # par(mfrow=c(3,1))
      # For plotting options, type '?plot.glmnet' in R console
       # plot(fit.ridge, xvar="lambda", main="Ridge")
       # plot(fit.lasso, xvar="lambda", main="LASSO")
       # plot(fit.elnet, xvar="lambda", main="Elastic Net")

       # par(mfrow=c(3,1))
       # plot(fit0, main="Ridge")
       # plot(fit1, main="LASSO")
       # plot(fit5, main="Elastic Net")
        # yhat        
        yhat1 <- predict(fit1, s=fit1$lambda.min, newx=x.test)
        yhat0 <- predict(fit0, s=fit0$lambda.min, newx=x.test)
        yhat5 <- predict(fit5, s=fit5$lambda.min, newx=x.test)
        # MSE
        mse_lasso[i] <- mean((y.test - yhat1)^2)
        mse_ridge[i] <- mean((y.test - yhat0)^2)
        mse_elanet[i] <- mean((y.test - yhat5)^2)	
}
	mse_lasso_final <- sum(mse_lasso)
      mse_ridge_final <- sum(mse_ridge)
      mse_elanet_final<- sum(mse_elanet)
      mse_lasso_final 
      mse_ridge_final 
      mse_elanet_final









      



