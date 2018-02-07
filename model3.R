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
chosen_variables <- c("Radio.Spend","Shopper.Marketing.Spend","Digital.FSI.Spend","Newspaper.FSI..Spend","WMX.Spend")
variables_set <- as.data.frame(t(combn(chosen_variables,2)))
colnames(variables_set) <- c("Variable1","Variable2")


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

result_dataframe3 <- data.frame(matrix(numeric(),ncol = 8, nrow = nrow(variables_set)))
colnames(result_dataframe3) <- c("model","variable0","variable1","variable2","summary","AIC","AICC","BIC")
for(i in 1 : nrow(variables_set)){
    tryCatch({
        print(i)
        model_num <- paste("model1",sep = "") 
        
        variable1 = data_div_sd[[variables_set[i,1]]]
        variable2 = data_div_sd[[variables_set[i,2]]]
        
        #It seems the example only uses tv as an independant variable. There are many more. How many should we attempt to incorporate in our model?
        #read prof. instructions
        
        ## Step 2. Find good starting values. To this end, let's apply Global Optimization
        
        
        ## Step 2A. Specify the model function
        
        nlsmean=function(b, zPriceVar, var1,var2) { #added parameters to input any variable for advertisement type and price
            #change here
            mu <- b[1]/zPriceVar +
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
        #change here 
        out_nls2 = nlsLM(zsales ~ b1/zprice +
                             (b2*exp(-b3*variable1))/(b4+exp(-b3*variable1)) +
                             (b5*exp(-b6*variable2))/(b7+exp(-b6*variable2)), start = list(b1=bb[1],b2=bb[2],b3=bb[3],b4=bb[4],b5=bb[5],b6=bb[6],b7=bb[7]))
        summary(out_nls2)
        
        
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
        
        result_dataframe3[[i,"model"]] <- model_num
        result_dataframe3[[i,"variable0"]] <- "price"
        result_dataframe3[[i,"variable1"]] <- as.character(variables_set[[i,1]])
        result_dataframe3[[i,"variable2"]] <- as.character(variables_set[[i,2]])
        result_dataframe3[[i,"summary"]]  <- list(summary(out_nls2))
        result_dataframe3[[i,"AIC"]] <- aic
        result_dataframe3[[i,"AICC"]] <- aicc
        result_dataframe3[[i,"BIC"]] <- bic },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# x <- data_div_sd$Radio
# y <-data_div_sd$Shopper.Marketing.Spend
# price <- mean(data_div_sd$Price)
# z_func <- function(price,adv_x,adv_y,b){
#     sales = b[1] + b[2]*exp(-b[3]*price) + b[4]*adv_x^b[5] + b[6]*adv_y^b[7]}
# z <- z_func(price,x,y,est)
# 
# plot_ly(data_div_sd, x = ~ x, y = ~ y, z = ~ z, type = 'scatter3d',opacity = 1)


##  Deliverables for HW 3
# 
# 1. R function to produce "My Estimates, Std Errors, t-values" 
# 
# 2. Try various nonlinear functions in Class3.pptx and various variables in class_data.csv. Apply information criterion to retain the best model
#   
# 3. Maximize Profit wrt to price and TV ad spends. Compare Actual vs Optimal and make your recommendations to the management team

