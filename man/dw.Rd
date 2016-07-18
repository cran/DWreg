\name{dw}
\alias{rdw}
\alias{qdw}
\alias{pdw}
\alias{ddw}
\title{Discrete Weibull}
\description{Density, distribution function, quantile function and random generation for the discrete Weibull distribution with parameters q and beta.}
\usage{
ddw(x,q=exp(-1),beta=1)
pdw(x,q=exp(-1),beta=1)
qdw(p,q=exp(-1),beta=1)
rdw(n,q=exp(-1),beta=1)
}
\arguments{
  \item{x}{quantile}
  \item{p}{probability}
  \item{n}{number of observations}
  \item{q,beta}{Parameters of the distribution}
}
\details{
The discrete Weibulll distribution has density
\deqn{p(x,q,\beta) = q^{x^{\beta}}-q^{(x+1)^{\beta}}}
  for \eqn{x = 0, 1, 2, \ldots}.
If q or beta are not specified they assume the default values of exp(-1) and 1, respectively. In this case, DW corresponds to a geometric distribution with p=1-q.  
  }
\value{
ddw gives the density, pdw gives the distribution function, qdw gives the quantile function, and rdw generates random samples from a DW distribution with parameters q and beta.
}
\references{Nagakawa T, Osaki S. The discrete Weibull distribution. IEEE transactions on reliability 1975; R-24(5).}
\author{Veronica Vinciotti}
\examples{
x<-rdw(1000,q=0.9,beta=1.5)
hist(x)
plot(x,unlist(lapply(x,ddw,q=0.9,beta=1.5)),ylab="density")
plot(x,unlist(lapply(x,pdw,q=0.9,beta=1.5)),ylab="cdf")
}
\keyword{dw}

