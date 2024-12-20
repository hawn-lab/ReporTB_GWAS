# Function to calculate model results for glm
glm_summ <- function(fit, SNP){
  #calculate sigma
  k=length(fit$coefficients)-1     #calculate the number of model parameters - 1
  SSE=sum(fit$residuals**2)        #calculate sum of squared residuals
  n=length(fit$residuals)          #calculate total observations in dataset
  sigma = sqrt(SSE/(n-(1+k)))      #calculate the residual standard error
  
  #format result
  model_results <- data.frame(tidy(fit),
                              AIC = fit$aic,
                              sigma= sigma,
                              snpID = SNP) %>% 
    rename(std_error=std.error, pval=p.value) %>% 
    select(term, snpID, estimate, std_error, pval, sigma, AIC)
  
  return(model_results)
}
