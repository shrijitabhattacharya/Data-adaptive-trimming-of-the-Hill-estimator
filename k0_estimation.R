censor_hill_plot = function(X,k0=0,k=max(k0+1,length(X)-1)){
  X.sort=sort(X,decreasing = TRUE);
  idx = c((k0+1):k);
  
  val = log(X.sort[(k0+1):(k+1)]);
  n = length(val);
  val = val[-n] - val[-1];
  val_wrong = val;
  val = cumsum(val*idx)/c(1:(k-k0));
  val_wrong = cumsum(val_wrong*c(1:(k-k0)))/c(1:(k-k0));
  
  return(list("xi"=val,"k0"=k0,"k"=k,"sigma"=val/sqrt(c(1:(k-k0))),
              "xi_wrong"=val_wrong,"sigma_wrong"=val_wrong/sqrt(c(1:(k-k0)))));
}

trim_hill_compute=function(X,k0,k){
  X.sort=sort(X,decreasing = TRUE)
  val=((k0+1)/(k-k0))*log(X.sort[(k0+1)]/X.sort[(k+1)])+(1/(k-k0))*sum(log(X.sort[(k0+2):k]/X.sort[(k+1)]))
  return (val)
}

trim_hill_plot = function(X,k){
  k0 = c(0:(k-1));
  xi = unlist(lapply(k0,function(x)trim_hill_compute(X,x,k)));
  return(list("k0"=k0,"xi"=xi,"sigma"=xi/sqrt(k-k0)));
}

compute_ADAP=function(X,k,level=0.05,a=1.2){
  T=rep(0,(k-1))
  E=unlist(lapply(0:(k-1), function(k0){trim_hill_compute(X,k0,k)}))
  U=rep(0,(k-1))
  C=rep(0,(k-1))
  count=1:(k-1)
  T=((k-count)/(k-count+1))*(E[(count+1)]/E[count])
  U=T^(k-count)
  C=(2*abs(U-0.5))[seq(from=(k-1),to=1,by=-1)]
  w0=a^(1:(k-1))
  w=(abs(log(1-level))*w0)/sum(w0)
  idx=which((abs(log(C))/w)<1)
  if (length(idx)==0){
    K0=0
  }else{
    K0=k-min(idx)
  }
  return(list('k0'=K0,'alpha'=1/E[(K0+1)],'xi'=E[(K0+1)]))
}



