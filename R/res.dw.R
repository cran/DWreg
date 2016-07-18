res.dw<-function(obj, k=5)
{
  y<-obj$data$response
  x<-obj$data$covariates[,-1]
  fitted.q<-obj$fitted.q
  fitted.beta<-obj$fitted.beta
  
  if(is.null(obj$call$para.q1))
    obj$call$para.q1<-FALSE
  if(is.null(obj$call$para.q2))
    obj$call$para.q2<-FALSE
  
  n=length(y)
  r.obs<-obj$residuals
  
  res.simenv=matrix(0,nrow=n,ncol=k)
  
  for(j in 1:k)
  {
    y.simenv=rdweibull(n,q=fitted.q, beta=fitted.beta,zero=TRUE)
    data.simenv=data.frame(y=y.simenv,x=x)
    fit.simenv=dw.reg(y~.,data=data.simenv,para.q1 = obj$call$para.q1,para.q2 = obj$call$para.q2,para.beta =obj$call$para.beta)
    
    fitted.q.simenv<-fit.simenv$fitted.q
    fitted.beta.simenv<-fit.simenv$fitted.beta
    
    a=1-fitted.q.simenv^(y.simenv^fitted.beta.simenv)
    b=1-fitted.q.simenv^((y.simenv+1)^fitted.beta.simenv)
    
    u=runif(n,a,b)
    r<-qnorm(u)
  
    res.simenv[,j]=sort(r)
  }
  
  quan<-t(apply(res.simenv,1,function(x) quantile(x,c(.025, .975))))
  rg=range(quan[,1],quan[,2],r.obs)
  
  n.point=qnorm((1:n)/(n+1))
  
  plot(n.point,sort(r.obs),ylim=rg,type="l",col=1, xlab="Standard normal quantiles",ylab = "Randomised quantile residuals" ,
       cex.lab=0.7,cex.axis=0.6,cex=0.6, main="q-q plot with 95% simulated envelope")
  
  polygon(c(n.point,rev(n.point)),c(quan[,1],rev(quan[,2])),col="gray87",border = NA)
  
  lines(n.point,sort(r.obs))
  
}
