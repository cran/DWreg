dw.reg<-function(formula, data,tau=0.5,para.q1=NULL,para.q2=NULL,para.beta=NULL, ...)
{
loglik.q1<-function(par,x,y)
{
beta<-par[length(par)]
if(beta<=0)
	loglik<-NA
else
{
theta<-par[-length(par)]
logq<-x%*%theta-log(1+exp(x%*%theta))
if(any(exp(logq)==0) | any(exp(logq)==1))
	loglik<-NA
else
	loglik<-sum(log(exp(y^beta*logq)-exp((y+1)^beta*logq)))
}
return(loglik)
}

loglik.beta<-function(par,x,y)
{
q<-par[length(par)]
if(q<=0 | q >=1)
	loglik<-NA
else
{
theta<-par[-length(par)]
beta<-exp(x%*%theta)
loglik<-sum(log(exp(y^beta*log(q))-exp((y+1)^beta*log(q))))
}
return(loglik)
}


loglik.q1beta<-function(par,x,y)
{
theta.q<-par[1:ncol(x)]
theta.beta<-par[(ncol(x)+1) : length(par)]
beta<-exp(x%*%theta.beta)
logq<-x%*%theta.q-log(1+exp(x%*%theta.q))
if(any(exp(logq)==0) | any(exp(logq)==1))
	loglik<-NA
else
	loglik<-sum(log(exp(y^beta*logq)-exp((y+1)^beta*logq)))
return(loglik)
}

loglik.q2beta<-function(par,x,y)
{
  theta.q<-par[1:ncol(x)]
  theta.beta<-par[(ncol(x)+1) : length(par)]
  beta<-exp(x%*%theta.beta)
  logq<--exp(x%*%theta.q)
  if(any(exp(logq)==0) | any(exp(logq)==1))
    loglik<-NA
  else
    loglik<-sum(log(exp(y^beta*logq)-exp((y+1)^beta*logq)))
  return(loglik)
}


if (is.null(para.q1) && is.null(para.q2) && is.null(para.beta)) 
  stop("Specify your model via para.q and para.beta!")

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
y <- model.response(mf, "numeric")
covariates <- model.matrix(mt, mf, contrasts)
term.labels <- colnames(covariates)



data<-data.frame(covariates,y)

if(is.null(para.q1))
  para.q1<-FALSE
if(is.null(para.q2))
  para.q2<-FALSE
if(is.null(para.beta))
  para.beta<-FALSE

oldw <- getOption("warn")
options(warn = -1)
q<-try(estdweibull(y,method="ML",zero=TRUE),silent=TRUE)[1]
beta<-try(estdweibull(y,method="ML",zero=TRUE),silent=TRUE)[2]
if(is.na(q) & para.beta)
  q<-try(estdweibull(y,method="P",zero=TRUE),silent=TRUE)[1]
if(is.na(q) & para.beta)
  q<-try(estdweibull(y,method="M",zero=TRUE),silent=TRUE)[1]
if(is.na(q) & para.beta)
{
  set.seed(131)
  q<-runif(1)
}

if(is.na(beta) & (para.q1 | para.q2))
  beta<-try(estdweibull(y,method="P",zero=TRUE),silent=TRUE)[2]
if(is.na(beta) & (para.q1 | para.q2))
  beta<-try(estdweibull(y,method="M",zero=TRUE),silent=TRUE)[2]
if(is.na(beta) & (para.q1 | para.q2))
{
  set.seed(131)
  beta<-runif(1)
}
options(warn = oldw)

if(para.q1){
int.in<-log(q/(1-q))
par.in<-c(int.in,rep(0,ncol(covariates)-1),beta)
names(par.in)<-c(term.labels,"beta")
mle<-maxLik(loglik.q1,start=par.in,x=covariates,y=y,...)
est<-coef(mle)
beta<-est[length(est)]
theta<-est[-length(est)]
logq<-covariates%*%theta-log(1+exp(covariates%*%theta))
Fittedy.tau<-ceiling((log(1-tau)/logq)^(1/beta)-1)
Fitted.q<-exp(logq)
Fitted.beta<-beta
}

if(para.q2){
  mle <- survreg(Surv(ifelse(y>0,y,NA),y+1,type="interval2")~covariates-1,data=data,dist="weibull")
  p<-length(coef(mle))
  medtab <- summary(mle)$table
  rownames(medtab) <- c(colnames(covariates),"Log(scale)")
  colnames(medtab) <- c("Estimate","Std.Err.","t value","Pr(>t)")
  medtab[,3] <- medtab[,1]/medtab[,2]
  medtab[,4] <- 2*pt(-abs(medtab[,3]), df=df.residual(mle))
  
  
  esttab <- matrix(nrow=p+1,ncol=4)
  rownames(esttab) <- c(colnames(covariates),"beta")
  colnames(esttab) <- c("Estimate","Std.Err.","t value","Pr(>t)")
  esttab[1:p,1] <- theta <- -coef(mle)/mle$scale
  esttab[p+1,1] <- beta <- 1/mle$scale
  esttab[,2] <- c(sqrt((coef(mle)/mle$scale)^2*((diag(mle$var)[1:p]/(coef(mle))^2)+(mle$var[p+1,p+1]))), sqrt((1/(mle$scale)^2)*mle$var[p+1,p+1]))
  esttab[,3] <- esttab[,1]/esttab[,2]
  esttab[,4] <- 2*pt(-abs(esttab[,3]),df=df.residual(mle))
  
  logq <- -exp(covariates%*%theta)
  Fittedy.tau<-ceiling((log(1-tau)/logq)^(1/beta)-1)
  Fitted.q<-exp(logq)
  Fitted.beta<-beta
  
}


if(para.beta)
{
int.in<-log(beta)
par.in<-c(int.in,rep(0,ncol(covariates)-1),q)
names(par.in)<-c(term.labels,"q")
mle<-maxLik(loglik.beta,start=par.in,x=covariates,y=y,...)
est<-coef(mle)
q<-est[length(est)]
theta<-est[-length(est)]
logbeta<-covariates%*%theta
Fittedy.tau<-ceiling((log(1-tau)/log(q))^(1/exp(logbeta))-1)
Fitted.q<-q
Fitted.beta<-exp(logbeta)
}

if(para.beta & para.q1)
{
int.in.beta<-log(beta)
int.in.q<-log(q/(1-q))
par.in<-c(int.in.q,rep(0,ncol(covariates)-1),int.in.beta,rep(0,ncol(covariates)-1))
names(par.in)<-c(term.labels,term.labels)
mle<-maxLik(loglik.q1beta,start=par.in,x=covariates,y=y,...)
est<-coef(mle)
theta.q<-est[1: (ncol(covariates))]
theta.beta<-est[(ncol(covariates)+1) : length(est)]
logq<-covariates%*%theta.q-log(1+exp(covariates%*%theta.q))
logbeta<-covariates%*%theta.beta
Fittedy.tau<-ceiling((log(1-tau)/logq)^(1/exp(logbeta))-1)
Fitted.q<-exp(logq)
Fitted.beta<-exp(logbeta)

}


if(para.beta & para.q2)
{
  int.in.beta<-log(beta)
  int.in.q<-log(-log(q))
  par.in<-c(int.in.q,rep(0,ncol(covariates)-1),int.in.beta,rep(0,ncol(covariates)-1))
  names(par.in)<-c(term.labels,term.labels)
  mle<-maxLik(loglik.q2beta,start=par.in,x=covariates,y=y,...)
  est<-coef(mle)
  theta.q<-est[1: (ncol(covariates))]
  theta.beta<-est[(ncol(covariates)+1) : length(est)]
  logq<--exp(covariates%*%theta.q)
  logbeta<-covariates%*%theta.beta
  Fittedy.tau<-ceiling((log(1-tau)/logq)^(1/exp(logbeta))-1)
  Fitted.q<-exp(logq)
  Fitted.beta<-exp(logbeta)
}

a=1-(Fitted.q)^(y^Fitted.beta)
b=1-(Fitted.q)^((y+1)^Fitted.beta)
u=runif(length(y),a,b)
res<-qnorm(u)

fit <- list()
class(fit) = "dw.reg"
fit$call <- call
fit$data<-list(covariates=data.frame(covariates),response=y)
fit$coefficients <- coef(mle)
fit$loglik <- logLik(mle)
fit$fitted.values <- Fittedy.tau
fit$fitted.q <- Fitted.q
fit$fitted.beta <- Fitted.beta
fit$residuals<-res

if(para.q2 & !para.beta)
{
  fit$tTable<-esttab
  fit$survreg<-medtab
}
else
  fit$tTable <- summary(mle)
return(fit)
}
