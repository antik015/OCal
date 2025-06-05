\name{posterior_projection}
\alias{posterior_projection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Orthogonal calibration for univariate outcomes.
}
\description{
This function implements orthogonal calibration when univariate (q=1) field measurements are available.
}
\usage{
posterior_projection(y, x, phi_sq, psi, nu, nmcmc, burnin, 
sig = sigT, m, prior = "BART", proj_type = "infinite", 
computer_model_functions, nugget = 0, alpha = 0.95, 
MH_sd = 0.1, prior_sd)
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
  \item{nmcmc}{
Number of MCMC iterations to perform.
}
  \item{burnin}{
Number of burnin samples.
}
  \item{sig}{
Plugin estimate of the error standard deviation.
}
  \item{m}{
Number of qudrature points to use for evaluating the projection. Usually, m = 15 works well.
}
  \item{prior}{
Prior choice for the bias function. Allows two choices "GP" and "BART".
}
  \item{proj_type}{
Type of projection. Allows two choices "infinite" and "finite".
}
  \item{computer_model_functions}{
A list of computer model functions. The list should contain the function definition along with derivatives of the function with respect to the calibration parameters.
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{thetaout}{Posterior samples of the calibration parameter. A p by (nmcmc - burnin) matrix when there are p calibration parameters.}
\item{time}{Time elapsed.}
\item{MH_acceptance}{Acceptance ratio of the Metropolis-Hastings algorithm for the calibration parameter.}
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
n = 100
sigT = 0.2
real_type = 1
sim_data = gen_data(n, real_type, sigT)
y = sim_data$y
x = sim_data$x
computer_model_functions = list()
computer_model_functions[[1]] = sim_data$f
computer_model_functions[[2]] = sim_data$g
nmcmc = 5000
burnin = 1000
phi_sq = 1
nu = 5/2
psi = 1/2
m = 15
nugget = 0
alpha = 0.95
MH_sd = 0.2 
prior_sd = 10
res1 = posterior_projection(y, x, phi_sq, nu, psi, nmcmc, burnin, sigT, m, "BART", "infinite", computer_model_functions, nugget, alpha, MH_sd, prior_sd)
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
