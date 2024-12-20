# Function to calculate model results for glmmkin
glmmkin_summ <- function(fit, SNP){
  #get model parameters
  beta <- fit$coefficients
  nvar <- length(beta)
  nfrail <- nrow(fit$cov) - nvar
  se <- sqrt(diag(fit$cov)[nfrail + 1:nvar])
  p.fit <- signif(1 - pchisq((beta/se)^2, 1), 2)
  AIC = length(fit$residuals)*(log(2*pi)+1+log((sum(fit$residuals^2)/length(fit$residuals)))) +
    ((length(fit$coefficients)+1)*2)
  
  #calculate sigma
  k=length(fit$coefficients)-1     #calculate the number of model parameters - 1
  SSE=sum(fit$residuals**2)        #calculate sum of squared residuals
  n=length(fit$residuals)          #calculate total observations in dataset
  sigma = sqrt(SSE/(n-(1+k)))      #calculate the residual standard error
  
  model_results <- data.frame(
    snpID = SNP,
    estimate = beta,
    std_error = se,
    pval = p.fit,
    sigma = sigma,
    AIC = AIC) %>% 
    rownames_to_column("term") %>% 
    select(term, snpID, estimate, std_error, pval, sigma, AIC)
  return(model_results)
}
