\name{dw.parest}
\alias{dw.parest}
\title{Parameter estimation for discrete Weibull}
\description{Estimation of the parameters q and beta of a discrete Weibull distribution}
\usage{
dw.parest(data,method,method.opt)
}
\arguments{
  \item{data}{Vector of observations}
  \item{method}{Either "likelihood" or "proportion"}
  \item{method.opt}{Optimization criterion used in maxLik (default is "NR")}
}
\details{If method="likelihood", the parameters q and beta are estimated by maximum likelihood. 

If method="proportion", the method of Araujio Santos and Fraga Alves (2013) is used, based on count frequencies.
}
\value{The function returns the parameter estimates of q and beta.}
\references{Araujo Santos P, Fraga Alves M. Improved shape parameter estimation in a discrete Weibull model. Recent Developments in Modeling and Applications in Statistics . Studies in Theoretical and Applied Statistics. Springer-Verlag, 2013; 71-80.}
\author{Veronica Vinciotti}
\examples{
x<-rdw(1000,q=0.9,beta=1.5)
dw.parest(x) #maximum likelihood estimates
dw.parest(x,method="proportion") #proportion estimates
}
\keyword{dw.parest}

