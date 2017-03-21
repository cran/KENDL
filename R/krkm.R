
silverman<-function(x)
{
  abs(0.5*exp(-abs(x)/sqrt(2))*sin(abs(x)/sqrt(2)+pi/4))/1.140086
}




ker<-function(obs,bdl,lod)

{
  n=length(lod)
  h=1.059*sd(lod)*n^(-1/3)
  wei=rep(1,length(obs))
  censored<- !bdl
  obs   <-   pmax(obs,lod)
  ind<-order(obs)
  data=cbind(obs,censored)
  data=data[ind,]
  obs=obs[ind]
  we=wei[ind]
  wd=lod[ind]
  dn<-data[,2]


  if(length(unique(obs))==length(obs)){
  ttn=rep(1,length(obs))
  een=dn
  temp=t(kronecker(t(wd),rep(1,n))-wd)
  temp=silverman(  (temp/h) )
  temp1=diag(dn)
  temp2=upper.tri(matrix(1,n,n), diag =TRUE )+matrix(0,n,n)
  temp3=matrix(1,n,n)
  temp6=sapply(obs,FUN=function(x) obs>=x)+matrix(0,n,n)
  num=temp%*%temp1
  deno=temp%*%temp2
  pro=1-(num/deno)
  cond=t(apply(pro,1,function (x){
    out=rev(cumprod(rev(x[2:n])))
    out=c(out,1)
    return(out)
  }  ))
  prob=apply(cond,2,mean,na.rm=T)

  temp7=matrix(diag(cond),n,n,byrow=FALSE)
  temp8=lower.tri(matrix(1,n,n), diag =FALSE )+matrix(0,n,n)
  temp9=matrix(prob,n,n,byrow=TRUE)
  inter_var=temp8*temp7+temp2*cond
  temp10=matrix(diag(cond),n,n,byrow=FALSE)
  temp55=matrix(dn,n,n,byrow=FALSE)
  variance<-cond-temp9-cond*((temp55*temp6)/(temp10)+1-(1/(inter_var)))
  sd=sqrt(apply(variance^2,2,sum,na.rm=TRUE))/n
  lower.cl=pmax(rep(0,n),prob-1.96*sd)
  upper.cl=pmin(rep(1,n),prob+1.96*sd)
  out<-cbind(unique.obs=obs,cdf=prob,se=sd,lower=lower.cl,upper=upper.cl)

  }else{
    ############################## tie   #######################
    unique.obs=unique(obs)
    een=tn=en=ttn=rep(NA,length(unique.obs))
    for (i in 1:length(unique.obs))
    {
      ttn[i]=sum(obs==unique.obs[i])
      tn[i]=sum(  we[obs==unique.obs[i]] )
      indi=which(obs==unique.obs[i])
      en[i]=sum(dn[indi]*we[indi])
      een[i]=sum(dn[which(obs==unique.obs[i])])
    }

    record=rep(NA,length(unique.obs))
    num=deno=temp.num=temp.deno=matrix(0,length(obs),length(unique.obs))
    for(i in 1: length(obs)){
      if(obs[i]<=max(unique.obs)){
        record=min(which(unique.obs>=obs[i]))
        temp.deno[i,record:length(unique.obs)]=1
        if(dn[i]==1){
          temp.num[i,record]=1
        }
      }
    }
    dist=t(kronecker(t(wd),rep(1,n))-wd)/h
    temp=silverman(dist)
    temp.num=temp.num
    temp.deno=temp.deno
    for(i in 1:length(lod))
    {
      for(j in 1:length(unique.obs))
      {
        num[i,j]=sum(temp[i,]*temp.num[,j])
        deno[i,j]=sum(temp[i,]*temp.deno[,j])
      }
    }
    pro=1-(num/deno)
    cond=t(apply(pro,1,function (x){
      out=rev(cumprod(rev(x[2:length(unique.obs)])))
      out=c(out,1)
      return(out)
    }  ))
    prob=apply(cond,2,mean)
    ind_lod=match(unique.obs,obs)
    cond=cond[ind_lod,]
    n1=length(ind_lod)
    temp7=matrix(diag(cond),n1,n1,byrow=FALSE)
    temp8=lower.tri(matrix(1,n1,n1), diag =FALSE )+matrix(0,n1,n1)
    temp9=matrix(prob,n1,n1,byrow=TRUE)
    temp2=upper.tri(matrix(1,n1,n1), diag =TRUE )+matrix(0,n1,n1)
    inter_var=temp8*temp7+temp2*cond
    temp10=matrix(diag(cond),n1,n1,byrow=FALSE)
    temp55=matrix(een,n1,n1,byrow=FALSE)
    temp6=sapply(unique.obs,FUN=function(x) unique.obs>=x)+matrix(0,n1,n1)
    variance<-cond-temp9-cond*((temp55*temp6)/(temp10)+1-(1/(inter_var)))
    sd=sqrt(apply(variance^2,2,sum,na.rm=TRUE))/n
    lower.cl=pmax(rep(0,n1),prob-1.96*sd)
    upper.cl=pmin(rep(1,n1),prob+1.96*sd)
    out<-cbind(unique.obs=unique.obs,cdf=prob,se=sd,lower=lower.cl,upper=upper.cl)
  }
  return(out)

}


KRKM <- function(obs,bdl,lod,method='formula',b=1000) UseMethod("KRKM")


KRKM.default<-function(obs,bdl,lod,method='formula',b=1000)
{
  formula=ker(obs,bdl,lod)
  if(method=='formula'){
    out=formula
  }else if (method=='bootstrap'){
    temp=matrix(NA,dim(formula)[1],b)
    for(i in 1:b)
    {
      ind=sample(1:length(lod),replace=TRUE)
      obs.b=obs[ind]
      lod.b=lod[ind]
      bdl.b=bdl[ind]
      x=formula[,c(1,2)]
      y=ker(obs.b,bdl.b,lod.b)[,c(1,2)]
      temp[,i]=merge(x,y,by='obs',all.x=TRUE)[,3]
    }
    se=apply(temp,1,sd,na.rm=TRUE)
    lower=pmax(rep(0,length(se)),out[,2]-1.96*se)
    upper=pmin(rep(1,(length(se))),out[,2]+1.96*se)
    var=cbind(se,lower,upper)
    out=formula[,1:2]
    out=cbind(out,var)
  }

  output=list()
  output$unique.obs=out[,1]
  output$cdf=out[,2]
  output$se=out[,3]
  output$lower=out[,4]
  output$upper=out[,5]
  class(output) <- 'KRKM'
  output
}

plot.KRKM<-function(x,conf.int=TRUE,lty=1,col=1,lwd=1,xlim=NULL,ylim=NULL,log="x",xlab=NULL,ylab='CDF',...)
{

  options( warn = -1 )

    plot(x$unique.obs,x$cdf,type='s',col=col,lty=lty,log=log,ylab=ylab,xlab=xlab,lwd=lwd,xlim=xlim,ylim=ylim,...)
    if(conf.int){
      points(x$unique.obs,x$lower,type='s',col='black',lty=2,log="x",lwd=1)
      points(x$unique.obs,x$upper,type='s',col='black',lty=2,log="x",lwd=1)
    }
}
