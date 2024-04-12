\name{posterior_projection_q}
\alias{posterior_projection_q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Orthogonal calibration for multivariate outcomes.
}
\description{
This function implements orthogonal calibration when multivariate (q>1) field measurements are available.
}
\usage{
posterior_projection_q(Y, x, phi_sq, nmcmc, burnin, Sig = SigT, m, prior = "BART", computer_model_functions, nugget = 0, alpha = 0.95, MH_sd = 0.1, prior_sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A n by q matrix of field observations.
}
  \item{x}{
A matrix with n rows of locations where field observations are available. Typical applications will have locations in a two-dimensional space.
}
  \item{phi_sq}{
Scale parameter of Gaussian process covariance if prior on bias function is GP. Typically, set to 1.
}
  \item{nmcmc}{
Number of MCMC iterations to perform.
}
  \item{burnin}{
Number of burnin samples.
}
  \item{Sig}{
Plugin estimate of the error covariance matrix.
}
  \item{m}{
Number of qudrature points to use for evaluating the projection. Usually, m = 15 works well.
}
  \item{prior}{
Prior choice for the bias function. Allows two choices "GP" and "BART".
}
  \item{computer_model_functions}{
A list of q computer model functions. Each list element should be another list containing the function definition along with derivatives of the function with respect to the calibration parameters.
}
  \item{nugget}{
Nugget to use if prior is "GP".
}
  \item{alpha}{
For finite projection, the weight of the sample covariance matrix while computing a linear shrinkage estimator of the population covariance matrix. Default is 0.95.
}
  \item{MH_sd}{
Proposal standard deviation for Random-walk Metropolis-Hastings for updating the calibration parameter.
}
  \item{prior_sd}{
Standard deviation of the Gaussian prior on the calibration parameter. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{thetaout}{Posterior samples of the calibration parameter. A p by (nmcmc - burnin) matrix when there are p calibration parameters.}
\item{time}{Time elapsed.}
\item{MH_acceptance}{Acceptance ratio of the Metropolis-Hastings algorithm for the calibration parameter.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Chakraborty, Antik, Jonelle B. Walsh, Louis Strigari, Bani K. Mallick, and Anirban Bhattacharya. "Orthogonal calibration via posterior projections with applications to the Schwarzschild model." arXiv preprint arXiv:2404.03152 (2024).
}
\author{
Antik Chakraborty
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
