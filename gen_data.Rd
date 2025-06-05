\name{gen_data}
\alias{gen_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate data for computer model calibration
}
\description{
This function generates data according to 3 different models given below.

1. For Model 1 the data is generated according to
\deqn{y_F(x) = y_R(x) + \epsilon, \, y_R(x) = 4x + x\sin 5x, \, f(x, t) = tx, \, \epsilon \sim N(0, \sigma^2), x \in [0,1]}

2. For Model 2 the data is generated according to
\deqn{y_F(x) = y_R(x) + \epsilon, \, f(x, t) = 7\{\sin (2\pi t_1 - \pi)\}^2 + 2\{(2\pi t_2 - \pi)^2 \sin (2\pi x - \pi)\}}
\deqn{y_R(x) = f(x, (0.2, 0.3)) , \, \epsilon \sim N(0, \sigma^2), x \in [0,1].}


3. Model 3 refers to the Friedman function. Relevant quantities are
\deqn{y_F(x) = y_R(x) + \epsilon, f(x, t) = t_1\sin x_1 x_2 + t_2(x_3 - 0.5)^2 + 10 x_4 + 5x_5}
\deqn{y_R(x) = f(x, (10, 20)), \epsilon \sim N(0,\sigma^2), x\in [0,1]^5.}
}
\usage{
gen_data(n, real_type, sigT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Sample size of observed (field) data. A positive integer.
}
  \item{real_type}{
Type of model. Can be either 1, 2, or 3 according to the above four models.
}
  \item{sigT}{
Value of \eqn{\sigma^2}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list with elements y, x, f (computer model), and associated derivatives g (if there are two calibration parameters returns g1 and g2).
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
real_type = 1
sigT = 0.2
d = gen_data(n, real_type, sigT)


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
