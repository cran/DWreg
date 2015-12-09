\name{dw.meanvar}
\alias{dw.meanvar}
\title{Mean and Variance of Discrete Weibull}
\description{Mean and variance of a discrete Weibull distribution with parameters q and beta.}
\usage{
dw.meanvar(q,beta,M)
}
\arguments{
  \item{q,beta}{Parameters of the distribution}
  \item{M}{Maximum value of the summation. Default value is 1000.}
}
\details{
The mean and variance are computed using the following approximations:
\deqn{E(X)=\sum_{k=1}^{M} q^{k^{\beta}}}
\deqn{E(X^2)=\sum_{k=1}^{M} (2k-1)q^{k^{\beta}} = 2\sum_{k=1}^{M} kq^{k^{\beta}}-E(X)}
}
\value{The function returns the mean and variance of a DW distribution with parameters q and beta.}
\references{Khan M, Khalique A, Abouammoth A. On estimating parameters in a discrete Weibull distribution. IEEE transactions on Reliability 1989; 38(3):348-350.}
\author{Veronica Vinciotti}
\examples{
dw.meanvar(q=0.9,beta=1.5)
#compare with sample mean/variance from a random sample
x<-rdw(1000,q=0.9,beta=1.5)
mean(x)
var(x)
}
\keyword{dw.meanvar}

