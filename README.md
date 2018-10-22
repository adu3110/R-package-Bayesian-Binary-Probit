# R-package-Bayesian-Binary-Probit

###This R package implements the Bayesian Auxiliary variable Model for Binary Regression

##Reference : Holems, C., Held, L. Bayesian Auxiliary Variable Models for 

##Binary and Multinomial Regression

##You can install the package from github by using the devtools package. 

##If you want the tar.gz file, please write a mail to the author of the package <aditi.jec31@gmail.com>. 

##Usage Example ##
```
library(devtools)

install_github("adu3110/R-package-bayesian-binary-probit")

library(BayesianBinaryProbit)

set.seed(250)
```
# Create 1000 samples of two independent variables
```
N <- 1000

X1 <- gaussiansamplesbyCLT(N)

X2 <- gaussiansamplesbyCLT(N)
```
# Create the matrix of input variables
```
X <- matrix(c(rep(1, N), X1, X2), ncol = 3)
```
# True values of regression coeffiecients
```
true_coefficients <- c(-0.8, 1.5, 0.6)
```
# Obtain the vector with probabilities of success p using the probit link
```
p <- pnorm(X %*% true_coefficients)
```
# Generate binary observation data y
```
y <- rbinom(N, 1, p)
```
# Fit the MLE Generalized Linear Model with probit link
```
fit <- glm(y ~ X1 + X2, family = binomial(link = probit))

fit$coefficients
```
#Create input data frame for Bayesian Model 
```
input_frame <- data.frame(resp = y, ind_var1 = X1, ind_var2 = X2)

system.time(

  bayesian_fit <- bayesianbinaryprobit(resp ~ ind_var1 + ind_var2, 
  
                                       data = input_frame, covar_prior = diag(10, 3), 
                                       
                                       num_mcmc = 20000, burn_in = 5000, thinning = 10)
                                       
)

bayesian_fit



##################Other Functions in the package##

##Gaussian Samples ##

hist(gaussiansamplesbyCLT(num_samples = 400, num_uniform_samples = 100, 

                     mean_norm = 5, sd_norm = 3))

#One sided truncated gaussian distribution

hist(samplegaussianonesided(800, side = "left", num_uniform_samples = 50))

hist(samplegaussianonesided(800, cut_off = 0.5, mean_norm = 2, sd_norm = 1,

                            side = "right", num_uniform_samples = 50))
```
#Inverse of a matrix by gauss elimination
```
gausseliminationinverse(matrix(c(2,1,1,0), 2, 2))
```
#Cholesky Lower Decomposition
```
choleskylower(matrix(c(2,-1,0,-1,2,1,0,1,2), 3, 3))
```
