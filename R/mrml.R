# 'weight' argument is not supported, use the same sampling strategy as for randomForest
# 11/2/2017
# 08/15/2018: added comparators per "ML Comparators".
# 12/4/2018: pass model formulas for R, A, Y models. Also add RYZ model.
# 04/30/2019: changed from mr_v2.R for package causalML.

mrml <- function(data, ind=1:nrow(data), model_R, model_A, model_Y, 
	Y_type="binary", tol=5e-3) {
# This algorithm is for a binary outcome Y.
# data: a data set, A=treatment, Y=outcome, R=1 when Z is observed, otherwise 0, X
#	is observed covariates.  Users have to use A, Y, R variable names for
#	treatment, outcome and missing indicator in data. 
# ind: For analysis using data as is, use the default.
#	the index of rows for the analysis, if ind=1:nrow(data), causal effect
#	for data as is.  It can be indices from sample with replacement.
# Y_type: outcome Y is either "binary" or "continuous"
# tol: for truncation of missing prob or propensity score to (tol, 1-tol).
	n <- nrow(data)
	bdat <- data[ind,]
	R <- bdat$R
# fit P(R=1|A,X,Y)
	R.fit <- rpart::rpart(formula=model_R, data=bdat, method="class",
		control=rpart.control(cp=0,xval=0))
	bdat$p_R <- predict(R.fit, type="prob")[,2]
	bdat$pi_tol <- with(bdat,ifelse(p_R < tol, tol, p_R))
# fit P(A=1|R=1, X, Z) with weights 1/p_R, samling to reflect the weighting
        prob_wt <- with(bdat,1/pi_tol/sum(1/pi_tol))
        N_w <- round(sum(1/bdat$pi_tol))
        bdat_1 <- bdat[sample(1:n, N_w, replace=T, prob = prob_wt),] 

        A.fit <- rpart::rpart(formula=model_A, data=bdat_1, subset=R==1,
                 method="class", control=rpart.control(cp=0,xval=0))
        bdat$p_A <- predict(A.fit, newdata=bdat, type="prob")[,2]

# fit P(Y=1|R=1, A, X, Z) with weights 1/p_R, predict for bdat
	  method <- ifelse(Y_type == "binary", "class","anova")
        Y.fit <- rpart::rpart(formula=model_Y, data=bdat_1, subset=R==1,
                method=method, control=rpart.control(cp=0,xval=0))
	  rm(bdat_1)

        bdat0 <- bdat
	  bdat0$A <- 0   
	  if (Y_type == "binary") {
	  	bdat$p_Y0 <- predict(Y.fit, newdata=bdat0, type="prob")[,2]
        } else bdat$p_Y0 <- predict(Y.fit, newdata=bdat0, type="vector")
    	  rm(bdat0)
  
        bdat1 <- bdat
	  bdat1$A <- 1
	  if (Y_type == "binary") {
	  	bdat$p_Y1 <- predict(Y.fit, newdata=bdat1, type="prob")[,2]
	  } else bdat$p_Y1 <- predict(Y.fit, newdata=bdat1, type="vector")
	  rm(bdat1)

# fit D model: 
# truncate p_A to (tol, 1-tol)
	bdat$e_tol <- ifelse(bdat$p_A< tol, tol, ifelse(bdat$p_A>1-tol, 1-tol, bdat$p_A))
	bdat$omega_0 <- with(bdat, (1-A)*Y/(1-e_tol)-((1-A)/(1-e_tol)-1)*p_Y0)
	bdat$omega_1 <- with(bdat, A*Y/e_tol-(A/e_tol-1)*p_Y1)

	D0.fit <- rpart::rpart(formula=update(model_R, omega_0~.), data=bdat, subset=R==1,
		method="anova", control=rpart.control(cp=0,xval=0))
	d0 <- predict(D0.fit, newdata=bdat, type="vector")
	D1.fit <- rpart::rpart(formula=update(model_R, omega_1~.), data=bdat, subset=R==1,
		method="anova", control=rpart.control(cp=0,xval=0))
	d1 <- predict(D1.fit, newdata=bdat, type="vector")

# phi, the causal effect

	phi <- with(bdat, mean(R/pi_tol*(omega_1-omega_0)-(R/pi_tol-1)*(d1-d0)))
# comparators per Changyu
# R and A models only:
	phi_ra <-with(bdat, mean(R/pi_tol*(A/e_tol - (1-A)/(1-e_tol))*Y))	
# R and Y models only:
	phi_ry <-with(bdat, mean(R/pi_tol*(p_Y1-p_Y0)))	
# R, A and Y models:
	phi_ray <-with(bdat, mean(R/pi_tol*(omega_1-omega_0)))	
# R, Y and Z models: just regress the mean model b(A,X,Z) (instead of DR) on A, X, Y
	b0.fit <- rpart::rpart(formula=update(model_R, p_Y0~.), data=bdat, subset=R==1,
		method="anova", control=rpart.control(cp=0,xval=0))
	b0 <- predict(b0.fit, newdata=bdat, type="vector")
	b1.fit <- rpart::rpart(formula=update(model_R, p_Y1~.), data=bdat, subset=R==1,
		method="anova", control=rpart.control(cp=0,xval=0))
	b1 <- predict(b1.fit, newdata=bdat, type="vector")

	phi_ryz <-with(bdat, mean(b1-b0))	
	c(phi, phi_ra, phi_ry, phi_ray, phi_ryz)
}
