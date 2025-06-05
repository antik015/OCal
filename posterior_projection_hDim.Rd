\name{posterior_projection_hDim}
\alias{posterior_projection_hDim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Orthogonal calibration for univariate outcomes with high-dimensional input.
}
\description{
This function is almost the same as posterior_projection(). The key difference here is that necessary integrals to evaluate the projection is computed using Monte Carlo instead of quadrature.
}
\usage{
posterior_projection_hDim(y, x, phi_sq, psi, nu, m, nmcmc, burnin, 
sig = sigT, prior = "BART", proj_type = "infinite", 
computer_model_functions, nugget = 0, MH_sd = 0.1, prior_sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A n by 1 vector of field observations.
}
  \item{x}{
A matrix with n rows of locations where field observations are available. 
}
  \item{phi_sq}{
Variance parameter for Matern kernel. 
}
  \item{psi}{
Lengthscale parameter for Matern kernel.
}
  \item{nu}{
Smoothness parameter of Matern kernel.
}
  \item{m}{
Number of Monte Carlo samples. m = 200 is a reasonable choice.
}
  \item{nmcmc}{
Number of MCMC iterations to perform.
}
  \item{burnin}{
Number of burnin samples.
}
  \item{sig}{
Plugin estimate of the error standard deviation.
}
  \item{prior}{
Prior choice for the bias function. Allows two choices "GP", "BART" and "BMARS".
}
  \item{proj_type}{
Type of projection. Supports "infinite" projections only.
}
  \item{computer_model_functions}{
A list of computer model functions. The list should contain the function definition along with derivatives of the function with respect to the calibration parameters.
}
  \item{nugget}{
Nugget to use if prior is "GP".
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
\item{theta_tilde}{Point at which orthogonality is enforced.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
n = 100
sigT = 0.2
real_type = 3
sim_data = gen_data(n, real_type, sigT)
y = sim_data$y
x = sim_data$x
computer_model_functions = list()
computer_model_functions[[1]] = sim_data$f
computer_model_functions[[2]] = sim_data$g1
computer_model_functions[[3]] = sim_data$g2
nmcmc = 5000
burnin = 1000
sig = sigT
prior_sd = 10
MH_sd = 0.05
prior = "BART"
proj_type = "infinite"
nugget = 0
phi_sq = 1
nu = 5/2
psi = 1/2
m = 200
res1 = posterior_projection_hDim(y, x, phi_sq, psi, nu, m, nmcmc, burnin, 
sigT, prior, proj_type, computer_model_functions, nugget,  
MH_sd, prior_sd)
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
