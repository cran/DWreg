\name{res.dw}
\alias{res.dw}
\title{DW regression: Diagnostics}
\description{Quantile-Quantile plot of the randomised quantile residuals of a DW regression fitted model with 95\% simulated envelope.}
\usage{res.dw(obj,k)
}
\arguments{
  \item{obj}{An object of class "dw.reg": the output of the dw.reg function.}
  \item{k}{The number of iterations for the simulated envelope.}
}
\details{Diagnostic check for a DW regression model. The randomised quantile residuals should follow a standard normal distribution.
}

\value{A q-q plot of the residuals with 95\% simulated envelope}
\references{Kalktawi, Vinciotti and Yu (2016) A simple and adaptive dispersion regression model for count data.}
\author{Veronica Vinciotti, Hadeel Kalktawi}
\examples{

#simulated example (para.q2=TRUE, beta constant)
theta0 <- -2
theta1 <- -0.5
beta<-0.5
n<-500
x <- runif(n=n, min=0, max=1.5)
logq<--exp(theta0  + theta1 * x)		
y<-unlist(lapply(logq,function(x,beta) rdw(1,q=exp(x),beta),beta=beta)) 
data.sim<-data.frame(x,y) #simulated data
fit<-dw.reg(y~x,data=data.sim,para.q2=TRUE)
res.dw(fit,k=5)
ks.test(fit$residuals,"pnorm")

#real example
library(Ecdat)
data(StrikeNb)
fit<-dw.reg(strikes~output,data=StrikeNb,para.q2=TRUE)
res.dw(fit,k=5)
ks.test(fit$residuals,"pnorm")
}
\keyword{res.dw}

