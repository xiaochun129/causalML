\name{mrml}
\alias{mrml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimating the average treatment effect in data subject to missing using tree models
}
\description{
mrml is used to compute the average treatment effect (ATE). It produces five ATE estimates which are based on five sets of model specifications:  phi based on (R, A, Y, Z), phi_ra based on (R, A), phi_ry based on (R, Y), phi_ray based on (R, A, Y, Z) and phi_ryz based on (R, Y, Z), where A is the binary treatment indicator (0/1), Y can be either a binary outcome (0/1) or a continuous outcome, Z is a covariate vector that is subject to missing, R=1 when Z is observed, otherwise 0, X is a fully observed covariate vector.
}
\usage{
mrml(data, ind=1:nrow(data), model_R, model_A, model_Y, Y_type="binary", tol=5e-3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame, containing the variables in the model.}
  \item{ind}{a vector of row indices. For analysis using data as is, use the default. It can used to either specify a subset of observations to be used in the analysis or a bootstrap sample.}
  \item{model_R, model_A, model_Y}{objects of class "formula": symbolic description of the models to be fitted for the indicator R for when Z is observed, for the treatment A, and for the outcome Y. The details of model specification are given under ‘Details’.}
  \item{Y_type}{a string, to specify that outcome Y is either "binary" or "continuous".}
  \item{tol}{a parameter for truncating P(R=1|A, X, Y) and the propensity score P(A=1|X, Z) to (tol, 1-tol) to satisfy the positivity assumption.}
}
\details{
Users have to use A, Y, R variable names for treatment, outcome and missing indicator in data. Fully observed covariates and covariates that are partially missing are user named. For example, the formulae for the R, A and Y models can be specified model_R <- R ~ A+age+bun+afib_01+creatinine+Y, model_A <- A ~ age+bun+afib_01+creatinine+queryEF, model_Y <- Y ~ A+age+bun+afib_01+creatinine+queryEF.
}
\value{mrml returns a vector of five estimates of ATE, phi_rayz, phi_ra, phi_ry, phi_ray, phi_ryz, based on model specifications of (R,A,Y,Z), (R,A), (R,Y), (R,A,Y) and (R,Y,Z) respectively.
}
\references{
Xiaochun Li, Zuoyi Zhang, Jinghua He, Shaun Grannis and Changyu Shen.  A Machine Learning Approach to Causal Inference in the Presence of Missing Data.
}
\author{Xiaochun Li}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(rpart)
# One MC data set with sample size 5000 for setting 2 in the paper.

data(eg_5000)

# The true ATE for the simulated data of setting 2

truth <- -0.040823585

# specify missing indicator, assignment, outcome models

model_R <- R ~ A+X1+X2+x3+X5+Y
model_A <- A ~ X1+X2+x3+X5+Z
model_Y <- Y ~ A+X1+X2+x3+X5+Z

# analyze the data set, no bagging (NOT recommended)

set.seed(123)
est <- mrml(eg_5000, ind=1:nrow(eg_5000), model_R, model_A, model_Y, Y_type="binary")
names(est) <-  c("phi_rayz", "phi_ra", "phi_ry", "phi_ray", "phi_ryz")
# bias
est - truth

# bagged estimate based on B bootstrap samples (Recommended)

set.seed(123)

# specify the number of bootstraps for the bagged estimate of ATE:

B <- 100  # we set to 100 in this example for demonstrative purpose; in paper we set to 1000.

nobs <- nrow(eg_5000)
phi_b <- replicate(B, mrml(eg_5000, sample(nobs, replace=TRUE), model_R, model_A, model_Y,Y_type="binary"))
# phi_b: 5 by B matrix
phi_bag <- rowMeans(phi_b, na.rm=T)
names(phi_bag) <- c("phi_rayz", "phi_ra", "phi_ry", "phi_ray", "phi_ryz")
phi_bag

# bias
phi_bag-truth


\dontrun{
# SE by bootstraping the bootstrap
library(doParallel)
library(foreach)
registerDoParallel(cores=50)

# set the number of bootstrap for estimating SE
BB <- 500

set.seed(12102018)

ptm <- proc.time()
res <- foreach(sim = 1:BB) \%dopar\% {
        dat_s <- eg_5000[sample(nobs,replace=TRUE),]
	phi_b <- replicate(B, mrml(dat_s, sample(nobs, replace=TRUE), model_R, model_A, model_Y))
	# phi_b: 5 by B matrix
	phi_bag <- rowMeans(phi_b, na.rm=T)
	names(phi_bag) <- c("phi_rayz", "phi_ra", "phi_ry", "phi_ray", "phi_ryz")
	phi_bag
}
proc.time()-ptm

res_arry <- sapply(res, I, simplify="array")
ave <- apply(res_arry, 1, mean, na.rm=T)
sds <- apply(res_arry, 1, sd, na.rm=T)
rmse <- sqrt((ave-truth)^2+sds^2)

# results of bias, SE and root-mean squared error
print(cbind(bias=ave-truth, sds=sds, rmse=rmse))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
