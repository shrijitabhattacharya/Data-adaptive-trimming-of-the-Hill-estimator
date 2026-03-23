## Simulation from distributions 

### Simulations from varying heavy tailed distributions for $\xi>0$

simulate_pareto=function(n,alpha,x0){
  u=runif(n)
  X=x0/u^(1/alpha)
  return (X)
}

simulate_burr=function(n,eta,lam,tau){
  U=runif(n)
  X=(eta*(U^(-(1/lam))-1))^(-(1/tau))
  return (X)
}

simulate_frec=function(n,alpha){
  U=runif(n)
  X=(-log(U))^(-(1/alpha))
  return (X)
}

simulate_T=function(n,df.T){
  X=rep(0,n)
  count=1
  while (count<=n){
    u=rt(1,df.T)                                                                                                                                                           
    if (u>0){
      X[count]=u
      count=count+1
    }
  }
  return (X)
}


simulate_gpd=function(n,mu,sigma,xi){
  U=runif(n)
  X=(sigma/xi)*(U^(-xi)-1)
  return (X+mu)
}   
## Simulations from heavy tails for $\xi<0$

simulate_gpd=function(n,mu,sigma,xi){
  U=runif(n)
  X=(sigma/xi)*(U^(-xi)-1)
  return (X+mu)
}                                 
