library(deSolve)

###################################
# SV model: lifelong protection

SV_mod <- function(t, y, parms) {
  #define state of system
  S <- y[1]   # susceptible
  V <- y[2]   # protected
  
  #determine change of system
  with(as.list(parms), {
    nu_v <- nu
    mu <- nu*S + nu_v*V
    
    dS = +mu -lambda*S -nu*S
    dV = +lambda*S -nu_v*V
    
    res <- c(dS, dV)
    
    list(res)
  })
}
 
# times <- seq(0, 100, length=1001)
# parms <- c(mu=0.02, nu=0.02, lambda=0.0115)
# 
# inits <- c(S=100, V=0)
# 
# #run model
# out <- as.data.frame(ode(method="lsoda", inits, times, SP_mod, parms))
# 
# 
# #visualize change
# plot(out$time, out$V, ylim=c(0,100), ylab="portion of total population (%)", xlab="time (yrs)", type="l")
# lines(out$time, out$S, lty=2)
# 
# legend("topright", legend=c("Protected", "Susceptible"),
#        lty=1:2, cex=0.8, horiz=FALSE)
# 
# text(x=out$time[nrow(out)]*0.975, y=out$V[nrow(out)]+(0.05*100), labels=format(round(out$V[nrow(out)],1),nsmall=1), cex=0.8)
# text(x=out$time[nrow(out)]*0.975, y=out$S[nrow(out)]-(0.05*100), labels=format(round(out$S[nrow(out)],1),nsmall=1), cex=0.8)
# 
# title(main="SV model - lifelong protection")
