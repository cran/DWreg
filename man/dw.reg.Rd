\name{dw.reg}
\alias{dw.reg}
\title{DW quantile regression}
\description{Parametric quantile regression for discrete response data. The conditional distribution of the response given the predictors is assumed to be DW with parameters q and beta dependent on the predictors.}
\usage{dw.reg(formula, data,tau=0.5,para.q1=FALSE,para.q2=TRUE,para.beta=FALSE,...)
}
\arguments{
  \item{formula}{An object of class "formula": a symbolic description of the model to be fitted.}
  \item{data}{An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which dw.qr is called.}
  \item{tau}{Quantile value (default 0.5). This is used only to extract the conditional quantile from the fitted distribution.}
  \item{para.q1, para.q2}{logical flag. If TRUE, the model includes a dependency of q on the predictors, as explained below.}
  \item{para.beta}{logical flag. If TRUE, the model includes a dependency of beta on the predictors, as explained below.}
  \item{...}{Additional arguments to the maxLik function}
}
\details{
The conditional distribution of Y (response) given x (predictors) is assumed a DW(q(x),beta(x)). 

If para.q1=TRUE, 
\deqn{log(q/(1-q))=\theta_0+\theta_1 X_1+\ldots+\theta_pX_p.} 

If para.q2=TRUE, 
\deqn{log(-log(q))=\theta_0+\theta_1 X_1+\ldots+\theta_pX_p.} 

If para.q1=FALSE and para.q2=FALSE, then q(x) is constant.

If para.beta=TRUE,
\deqn{log(\beta)= \gamma_0+\gamma_1 X_1+\ldots+\gamma_pX_p.} Otherwise beta(x) is constant.

By default, para.q2 is set to TRUE and the rest to FALSE, as in the continuous Weibull regression model.
}

\value{ A list of class dw.qr containing the following components:

\item{call}{the matched call.}
\item{coefficients}{the theta and gamma coefficients as defined above.} 
\item{fitted.values}{fitted values (on the response scale) for the specified quantile tau.}
\item{tTable}{coefficients, standard errors, etc.} 
}
\references{Kalktawi, Vinciotti and Yu (2015) A simple and adaptive dispersion regression model for count data.}
\author{Veronica Vinciotti, Hadeel Kalktawi}
\examples{
#simulated example (para.q1=TRUE)
theta0 <- 2
theta1 <- 0.5
beta<-0.5
n<-500
x <- runif(n=n, min=0, max=1.5)
logq<-theta0 + theta1 * x - log(1+exp(theta0  + theta1 * x))		
y<-unlist(lapply(logq,function(x,beta) rdw(1,q=exp(x),beta),beta=beta)) 
data.sim<-data.frame(x,y) #simulated data
fit<-dw.reg(y~x,data=data.sim,para.beta=FALSE,para.q1=TRUE,para.q2=FALSE)
fit$tTable	

#simulated example (para.q2=TRUE)
theta0 <- -2
theta1 <- -0.5
beta<-0.5
n<-500
x <- runif(n=n, min=0, max=1.5)
logq<--exp(theta0  + theta1 * x)		
y<-unlist(lapply(logq,function(x,beta) rdw(1,q=exp(x),beta),beta=beta)) 
data.sim<-data.frame(x,y) #simulated data
fit<-dw.reg(y~x,data=data.sim,para.beta=FALSE,para.q1=FALSE,para.q2=TRUE)
fit$tTable	


#real example
library(Ecdat)
data(StrikeNb)
fit<-dw.reg(strikes~output,data=StrikeNb,para.beta=FALSE,para.q1=FALSE,para.q2=TRUE)
fit$tTable
}
\keyword{dw.reg}

