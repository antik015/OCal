\name{OCal-package}
\alias{OCal-package}
\alias{OCal}
\docType{package}
\title{
\packageTitle{OCal}
}
\description{
\packageDescription{OCal}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{OCal}
\packageIndices{OCal}
Given field observations and a computer code (explicitly known), statistical calibration estimates the best parameter combination of the computer model for approximating the physical process that generated the data. This package develops an orthogonal calibration method where idenfiability of the calibration parameter is ensured by a linear constriant on the so-called bias function. Unlike existing methods, the method makes the choice of the bias function free by using projected posteriors. For more details, refer to "Orthogonal calibration via posterior projections with applications to the Schwarzschild model" by Chakraborty et.al. (2024), which can be found at  https://arxiv.org/abs/2404.03152.
}
\author{
\packageAuthor{OCal}

Maintainer: \packageMaintainer{OCal}
}
\references{
Chakraborty, Antik, Jonelle B. Walsh, Louis Strigari, Bani K. Mallick, and Anirban Bhattacharya. "Orthogonal calibration via posterior projections with applications to the Schwarzschild model." arXiv preprint arXiv:2404.03152 (2024).
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
f1 = function(t, x)
{
  return(t*x)
}
y1_R = function(x)
{
  r = 4*x + x*sin(5*x)
  return(r)
}
f2 = function(t, x)
{
  r = pnorm(t*(x - 0.5))
  return(r)
}

y2_R = function(x)
{
  r = exp(6*(x - 0.5))/(1 + exp(6*(x - 0.5)))
  return(r)
}

gen_data2 = function(n, SigT)
{
  ## Computer model
  f1 = function(t, x)
  {
    return(t*x)
  }
  
  # f = function(x, t)
  # {
  #   return(t*x)
  # }
  
  ## Computer model derivative
  g11 = function(t, x)
  {
    return(x)
    #return(x*sin(t*x))
  }
  
  ## Physical process
  y1_R = function(x)
  {
    r = 4*x + x*sin(5*x)
    return(r)
  }
  
  f2 = function(t, x)
  {
    r = pnorm(t*(x - 0.5))
    return(r)
  }
  
  g21 = function(t, x)
  {
    r = dnorm(t*(x-0.5))*(x - 0.5)
    return(r)
  }
  
  y2_R = function(x)
  {
    r = exp(6*(x - 0.5))/(1 + exp(6*(x - 0.5)))
    return(r)
  }
  
  x = runif(n)
  M = matrix(0, n, 2)
  M[,1] = y1_R(x)
  M[,2] = y2_R(x)
  y = M + rmvnorm(n, sigma = SigT)
  computer_model = list()
  computer_model[[1]] = list("f1" = f1, "g11" = g11)
  computer_model[[2]] = list("f1" = f2, "g11" = g21)
  
  result = list("y" = y, "x" = x, "SigT" = SigT, "computer_model" = computer_model)
  return(result)
}

n = 100
s0 = 0.2
SigT = s0^2*matrix(c(1, 0.5, 0.5, 1), 2, 2)
d = gen_data2(n, SigT)

nmcmc = 5000
burnin = 1000
m = 15
phi_sq = 1
prior_sd = 10
MH_sd = 0.1

computer_model1 = d$computer_model[[1]]
computer_model2 = d$computer_model[[2]]


## Estimate sigma from both model ##
bart_fit1 = wbart(d$x, d$y[,1])
e1 = d$y[,1] - bart_fit1$yhat.train.mean

bart_fit2 = wbart(d$x, d$y[,2])
e2 = d$y[,2] - bart_fit2$yhat.train.mean
S = cov(cbind(e1, e2))

res_bart = posterior_projection_q(d$y, d$x, phi_sq, nmcmc, burnin, Sig = S, m, prior = "BART", 
                                computer_model_functions = d$computer_model, MH_sd = MH_sd, prior_sd = prior_sd)
}
