#setwd("C:/Users/sprot/Dropbox (School of Management)/UCD/442_stats/Class 3")
library(tidyverse)
library(numDeriv)
library(GenSA)
library(nlstools)
library(nls2)
library(corrplot)
library(plotly)


data=read.csv("class_data.csv", header=T)    # read csv file and label the data as "data"

corrplot(cor(data[,-1]), method = "color")

time=data[,1]
sales=data[,2]
price=data[,3]
tv=data[,4]
print=data[,5]
radio=data[,6]
psearch=data[,7]
psocial=data[,8]
shopper=data[,9]
ntdola=data[,10]
ntdolv=data[,11]
tdola=data[,12]
tdolv=data[,13]
email=data[,14]
linqia=data[,15]
linqia_shopper=data[,16]
dfsi=data[,17]
newsfsi=data[,18]
wmx=data[,19]

##### Step-by-Step Nonlinear Regression  #####
myDivSd = function(vecX){
    divBySd = vecX/sd(vecX)
}
zsales=sales/sd(sales)
zprice=price/sd(price)
ztv=tv/sd(tv)
zprint = print/sd(print)
data_div_sd = data %>%
    mutate_all(funs(myDivSd(.))) #does the division by sd for all columns
## Step 1. Scale variables. Most important step. Scale them such that parameter estimates acquire similar magnitudes

variable1 = data_div_sd[["Shopper.Marketing.Spend"]]
variable2 = data_div_sd[["WMX.Spend"]]

#It seems the example only uses tv as an independant variable. There are many more. How many should we attempt to incorporate in our model?
#read prof. instructions

## Step 2. Find good starting values. To this end, let's apply Global Optimization


## Step 2A. Specify the model function

nlsmean=function(b, zPriceVar, var1,var2) { #added parameters to input any variable for advertisement type and price
    #change here
    mu <- b[1] + b[2]*exp(-b[3]*zPriceVar) +
        (b[2]*exp(-b[3]*var1))/(b[4]+exp(-b[3]*var1)) +
        (b[5]*exp(-b[6]*var2))/(b[7]+exp(-b[6]*var2))		# mu denotes the mean function
    return(mu)
}

# slide8Model1 = function(b, var1,var2){
#     sales2 = b[1]+b[2]*exp(-b[3]*zprice) + b[4]*var1^b[5] + b[6]*var2^b[7]
#     return(sales2)
# }

## Step 2B. Specify SSE function

sse=function(b) {
    
    sse=sum((zsales-nlsmean(b, data_div_sd$Price , variable1,variable2))^2)
}

## Step 2C. Find good starting values using global minimizer like Simualted Annealing (or Genetic Algorithm or Particle Swarm)

par=c(1,1,1,1,1,1,1)								  # starting values for the global optimizer
lower=c(-15,-15,-3,-15,-3,-15,-3)  				# lower bounds on parameter values
upper=c(15,15,3,15,3,15,3)  					# upper bounds on parameter values
out=GenSA(par, sse, lower, upper)			# GenSA = Generalized Simulated Annealing bring us near the neighborhood of global solutions


# Step 3. Nonlinear Regression using built-in R (may not work)

bb=out$par										  # use GenSA solution as the starting values for nls()
#out_nls=nls(zsales~b1+b2/(zprice)^b3 + b4*(ztv)^b5,start=list(b1=bb[1],b2=bb[2],b3=bb[3],b4=bb[4],b5=bb[5]))


library(minpack.lm)
## Step 4. Nonlinear Regression from first principles (works often)


par=bb											   # use GenSA solution as the starting values for nls()
fit=optim(par, sse, method = "BFGS", hessian = T) 	# optim = solver that minimizes a specified fcn (here sse)
est=fit$par									# final estimates from optim solver

## Inference: Used for Confidence Intervals

nn=ncol(t(zsales))								# sample size
pp=ncol(t(est))								# number of parameters
yhat=nlsmean(est, zPriceVar = zprice, var1 = variable1, var2 = variable2)								# forecast zsales
err=zsales-yhat								# residuals
sig2=sum(err^2) /(nn-pp)						# sigma^2 of error term
jmat=jacobian(nlsmean,x = est, zPriceVar = zprice, var1 = variable1, var2 = variable2)						# jmat = gradient of the mean function at the estimated parameters
jmat
varp=sig2*solve(t(jmat) %*% jmat+0.1*diag(7)) 		# variance-covariance matrix of parameters. I added small ridge regularization	to ensure inverse
se=sqrt(diag(varp))							# diag elements are variances, sqrt makes them std dev
tvals=est/se									# tvals for inference: abs(tval) > 1. 65 => 90% confidence level; abs(tval) > 1.96 => 95% CI

## Information Criteria: Used for Model Selection

aic = nn*log(sig2) + 2*pp						# Use when nn is large (i.e., pp/nn < 5%)
aicc = nn*log(sig2) + nn*(nn+pp)/(nn-pp-2)			# Use when pp/nn > 5% => invented by Prof. Tsai at GSM and Hurvich at NYU
bic = nn*log(sig2) + pp*log(nn)					# Use when nn is large (i.e., pp/nn < 5%)
model_parameter <- data.frame(aic, aicc, bic)

est <- c(est)
tvals <- c(tvals)
se <- c(se)
df_result <- data.frame(est, tvals, se)

##############
#### HW 4 ####
##############
b1 = 0
# b1 = 1.634239 < 1.65 --> 0 --> 95%
b2 = df_result$est[2]
b3 = df_result$est[3]
b4 = df_result$est[4]
b5 = df_result$est[5]
b6 = df_result$est[6]
b7 = df_result$est[7]
c <- 1 # variable cost
profit <- function(x) {
  sales = b1 + b2*exp(-b3*x[1]) +
    (b2*exp(-b3*x[2]))/(b4+exp(-b3*x[2])) +
    (b5*exp(-b6*x[3]))/(b7+exp(-b6*x[3]))
  profit = (x[1]-c)*sales - x[2] - x[3] # x[1] is zprice, x[2] and x[3] are the fixed costs
  return(profit)
}
x0 <- c(10,10,10)
out = optim(x0, profit, method = "L-BFGS-B", hessian = TRUE) 
print("Max Profit, Optimal Price, Optimal Spend"); cbind(out$value,out$par[1], out$par[2], out$par[3])
