library(deSolve)

###################################
# SPR model: protection is lost over time

SPR_mod <- function(t, y, parms) {
  #define state of system
  S <- y[1]                       # susceptible
  V_0 <- y[2]                     # protected
  eval(parse(text = paste0("V_",parms["iN"]," <- y[",parms["iN"]+2,"]"))) 
  R <- y[parms["iN"]+3]           # resusceptible
  
  V <- sum(y[2:(parms["iN"]+2)])                        # total number of protected
  SR <- sum(y[c(1,parms["iN"]+3)])                      # total number of unprotected
  
  #determine change of system
  with(as.list(parms), {
    nu_v <- nu*f_nu
    mu <- nu*SR + nu_v*V
    
    dS = +mu -nu*S -lambda*S
    dV_0 = +lambda*S -gamma*V_0 -nu_v*V_0 
    res <- c(dS,dV_0)
    
    if(iN > 0){
      Phi <- diag(-gamma-nu_v,iN)
      OffDiag <- row(Phi) - col(Phi)
      Phi[OffDiag == 1] <- gamma

      dV_iN <- as.numeric(Phi %*% y[3:(iN+2)])
      res[3] = +gamma*V_0 + dV_iN[1]
      res[4:(iN+2)] = dV_iN[-1]
    }
    
    eval(parse(text = paste0("dR = -nu*R +gamma*V_",iN)))
    res[iN+3] <- dR
    
    list(res)
  })
}

# times <- seq(0, 100, length=101)
# parms <- c(mu=0.02, nu=0.02, lambda=0.05, f_nu=0.5)
# 
# #loss of protection
# m_dur <- 30 #mean time until protection is lost (yrs)
# var_dur <- 1.5^2 #variance of time until protection is lost (yrs)
# parms["gamma"] <- m_dur/var_dur - parms["nu"]
# parms["iN"] <- round(var_dur*(parms["gamma"]+parms["nu"])^2 -1)
# message(paste("approximated time until loss of protection is", round((parms["iN"]+1)/(parms["gamma"]+parms["nu"]),1),"years"))
# 
# inits <- c(S=100, V_0=0)
# if(parms["iN"] > 0){
#   for(i in 1:parms["iN"]) inits[paste0("I_",i)] <- 0 #add initial values of protected in transit
# }
# inits["R"] <- 0
# 
# 
# #run model
# system.time(out <- as.data.frame(ode(method="lsoda", inits, times, SPR_mod, parms)))
# 
# #visualize change
# V <- apply(out[3:(dim(out)[2]-1)],1,sum) #sum of all protected in transit
# plot(out$time, V, ylim=c(0,100), ylab="portion of total population (%)", xlab="time (yrs)", type="l")
# lines(out$time, out$S, lty=2)
# lines(out$time, out$R, lty=3)
# 
# legend("topright", legend=c("Protected", "Susceptible","Resusceptible"),
#        lty=1:3, cex=0.8, horiz=FALSE)
# 
# text(x=out$time[nrow(out)]*0.975, y=V[length(V)]+(0.05*100), labels=format(round(V[length(V)],1),nsmall=1), cex=0.8)
# text(x=out$time[nrow(out)]*0.975, y=out$S[nrow(out)]-(0.05*100), labels=format(round(out$S[nrow(out)],1),nsmall=1), cex=0.8)
# text(x=out$time[nrow(out)]*0.975, y=out$R[nrow(out)]+(0.05*100), labels=format(round(out$R[nrow(out)],1),nsmall=1), cex=0.8)
# 
# title(main=paste("SPR model - protection is lost over time\ntime until loss:",m_dur,"±",sqrt(var_dur),"years"))
# 
