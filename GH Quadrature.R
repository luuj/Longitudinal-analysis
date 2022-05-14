library(fastGHQuad)

#f: Function to integrate with respect to first (scalar) argument; this does not include the weight function exp(-x^2)
#rule: Gauss-Hermite quadrature rule to use, as produced by gaussHermiteData

labFunction <- function(beta0_star, beta1_star, sigma, mu, nnodes){
  rule <- gaussHermiteData(nnodes)
  
  expit <- function(x){plogis(x)}
  logit <- function(x){qlogis(x)}
  
  f <- function(gamma, z){
    pi^(-0.5) * expit(beta0_star + beta1_star*z + sqrt(2*sigma^2) * gamma + mu)
  }
  
  browser()
  
  beta1 <- logit(ghQuad(f, rule, z=1)) - logit(ghQuad(f, rule, z=0))
  
  return(exp(beta1))
}

labFunction(-2, 0.4, sqrt(2), 0, 10)

