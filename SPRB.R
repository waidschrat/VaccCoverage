library(deSolve)

###################################
# SPRB model: protection is lost over time
# but may be refreshed by a booster

SPRB_mod <- function(t, y, parms) {
  #define state of system
  S <- y[1]                       # susceptible
  V_0 <- y[2]                     # protected
  if(parms["iN"] > 0 & parms["iN"] >= parms["iB"]){ #time until boosting < time until protection loss
    eval(parse(text = paste0("V_",parms["iN"]," <- y[",parms["iN"]+2,"]")))
    B <- y[parms["iN"]+3]         # boosted
    R <- y[parms["iN"]+4]         # resusceptible
    
    SR <- sum(y[c(1,(parms["iN"]+3):(parms["iN"]+4))])  # total number of unprotected
    
  }else if(parms["iN"] > 0 & parms["iN"] < parms["iB"]){ #time until boosting > time until protection loss
    eval(parse(text = paste0("V_",parms["iB"]," <- y[",parms["iB"]+2,"]")))
    B <- y[parms["iB"]+3]         # boosted
    R <- y[parms["iB"]+4]         # resusceptible
    
    SR <- sum(y[c(1,(parms["iN"]+3):(parms["iB"]+4))])  # total number of unprotected
  }
  V <- sum(y[2:(parms["iN"]+2)])                        # total number of protected
  
  #determine change of system
  with(as.list(parms), {
    nu_v <- nu*f_nu
    mu <- nu*SR + nu_v*V
    
    dS = +mu -nu*S -lambda*S
    dV_0 = +lambda*S -gamma*V_0 -nu_v*V_0 +100*B
    res <- c(dS,dV_0)
    
    if(iN > 0 & iN >= iB){ #time until boosting < time until protection loss
      Phi <- diag(-gamma-nu_v,iN+1)
      OffDiag <- row(Phi) - col(Phi)
      Phi[OffDiag == 0][iB:iN] <- -gamma-nu_v-delta
      Phi[OffDiag == 1] <- gamma
      Phi[nrow(Phi),iB:iN] <- delta #booster input from transit compartments
      Phi[nrow(Phi),iN+1] <- -100 #booster output

      dV_iN <- as.numeric(Phi %*% y[3:(iN+3)])
      res[3] = +gamma*V_0 + dV_iN[1]
      res[4:(iN+2)] = dV_iN[c(-1,-(iN+1))]
      res[iN+3] = dV_iN[iN+1] +delta*R #booster compartment
      
      eval(parse(text = paste0("dR = -nu*R +gamma*V_",iN," -delta*R")))
      res[iN+4] <- dR
    }else if(iN > 0 & iN < iB){ #time until boosting > time until protection loss
      Phi <- diag(-gamma-nu,iB+1)
      OffDiag <- row(Phi) - col(Phi)
      Phi[OffDiag == 0][1:iN] <- -gamma-nu_v
      Phi[OffDiag == 0][iB] <- -gamma-nu-delta
      Phi[OffDiag == 1] <- gamma
      Phi[nrow(Phi),iB] <- delta #booster input from last transit compartment
      Phi[nrow(Phi),iB+1] <- -100 #booster output

      dV_iN <- as.numeric(Phi %*% y[3:(iB+3)])
      res[3] = +gamma*V_0 + dV_iN[1]
      res[4:(iB+2)] = dV_iN[c(-1,-(iB+1))]
      res[iB+3] = dV_iN[iB+1] +delta*R #booster compartment
      
      eval(parse(text = paste0("dR = -nu*R +gamma*V_",iB," -delta*R")))
      res[iB+4] <- dR
    }
    
    list(res)
  })
}


# times <- seq(0, 100, length=101)
# parms <- c(mu=0.02, nu=0.02, lambda=0.05, delta=0.08, f_nu=0.5)
# 
# #loss of protection
# m_dur <- 5 #mean time until protection is lost (yrs)
# var_dur <- 2^2 #variance of time until protection is lost (yrs)
# parms["gamma"] <- m_dur/var_dur - parms["nu"]
# parms["iN"] <- round(var_dur*(parms["gamma"]+parms["nu"])^2 -1)
# message(paste("approximated time until loss of protection is", round((parms["iN"]+1)/(parms["gamma"]+parms["nu"]),1),"years"))
# 
# #boosting interval
# boost <- 30 #time interval for vaccination boosting (yrs)
# parms["iB"] <- round((var_dur*(parms["gamma"]+parms["nu"])^2 -1)*(boost/m_dur))
# message(paste("approximated start time of boosting is", round((parms["iB"]+1)/(parms["gamma"]+parms["nu"]),1),"years"))
# 
# 
# inits <- c(S=100, V_0=0)
# if(parms["iN"] > 0 & parms["iN"] >= parms["iB"]){
#   for(i in 1:parms["iN"]) inits[paste0("I_",i)] <- 0 #add initial values of protected in transit
# }else if(parms["iN"] > 0 & parms["iN"] < parms["iB"]){
#   for(i in 1:parms["iB"]) inits[paste0("I_",i)] <- 0 #add initial values of protected in transit
# }
# inits["B"] <- 0
# inits["R"] <- 0
# 
# #run model
# system.time(out <- as.data.frame(ode(method="lsoda", inits, times, SPRB_mod, parms)))
# 
# #visualize change
# if(parms["iN"] >= parms["iB"]){
#   V <- apply(out[3:(dim(out)[2]-2)],1,sum) #sum of all protected in transit
#   R <- out$R #all resusceptible
# }else if(parms["iN"] < parms["iB"]){
#   V <- apply(out[3:(dim(out)[2]-2-parms["iB"]+parms["iN"])],1,sum) #sum of all protected in transit
#   R <- apply(out[(dim(out)[2]-1-parms["iB"]+parms["iN"]):(dim(out)[2]-2)],1,sum) + out$R #all resusceptible
# }
# 
# # dat <- data.frame(matrix(c(parms,m_dur,boost), nrow=1, dimnames = list(NULL,c(names(parms),"m_dur","boost"))))
# # plotSV(out,dat,1)
# 
# plot(out$time, V, ylim=c(0,100), ylab="portion of total population (%)", xlab="time (yrs)", type="l")
# lines(out$time, out$S, lty=2)
# lines(out$time, R, lty=3)
# 
# legend("topright", legend=c("Vaccinated", "Susceptible","Resusceptible"),
#        lty=1:3, cex=0.8, horiz=FALSE)
# 
# text(x=out$time[nrow(out)]*0.975, y=V[length(V)]+(0.05*100), labels=format(round(V[length(V)],1),nsmall=1), cex=0.8)
# text(x=out$time[nrow(out)]*0.975, y=out$S[nrow(out)]-(0.05*100), labels=format(round(out$S[nrow(out)],1),nsmall=1), cex=0.8)
# text(x=out$time[nrow(out)]*0.975, y=R[length(R)]+(0.05*100), labels=format(round(R[length(R)],1),nsmall=1), cex=0.8)
# 
# title(main=paste("SPRB model - boosting possible after:",round((parms["iB"]+1)/(parms["gamma"]+parms["nu"]),1),"years\ntime until protection loss:",m_dur,"±",sqrt(var_dur),"years"))
