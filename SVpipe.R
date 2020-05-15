PruneOutput <- function(dat,out){
  for(i in 1:length(out)){
    S <- out[[i]]$S
    
    if(dat[i,"m_dur"] == Inf){
      V <- out[[i]]$V
      R <- rep(0, length(out[[i]]$time))
    }else if(dat[i,"m_dur"] != Inf & dat[i,"boost"] == Inf){
      V <- apply(out[[i]][3:(ncol(out[[i]])-1)],1,sum) #sum of all protected in transit
      R <- out[[i]]$R
    }else{
      if(dat[i,"iN"] >= dat[i,"iB"]){
        V <- apply(out[[i]][3:(ncol(out[[i]])-2)],1,sum) #sum of all protected in transit
        R <- out[[i]]$R #all resusceptible
      }else if(dat[i,"iN"] < dat[i,"iB"]){
        V <- apply(out[[i]][3:(ncol(out[[i]])-2-dat[i,"iB"]+dat[i,"iN"])],1,sum) #sum of all protected in transit
        R <- apply(out[[i]][(ncol(out[[i]])-1-dat[i,"iB"]+dat[i,"iN"]):(ncol(out[[i]])-2)],1,sum) + out[[i]]$R #all resusceptible
      }
    }
    
    out[[i]] <- data.frame("time"=out[[i]]$time, "S"=S, "V"=V, "R"=R)
  }
  return(out)
}


SVpipe <- function(times, parms, init=NULL, method="lsoda", type="SP"){
  if(type == "SP"){
    if(parms["m_dur"] == Inf){
      assign("mod", SP_mod)
      inits <- c(S=100, V=0)
      if(!is.null(init)) inits[c("S","V")] <- init[c("S","V")]
      
    }else if(parms["m_dur"] != Inf & parms["boost"] == Inf){
      assign("mod", SPR_mod)
      inits <- c(S=100, V_0=0)
      if(!is.null(init)) inits[c("S","V_0")] <- init[c("S","V")]
      if(parms["iN"] > 0){
        if(!is.null(init)){
          for(i in 0:parms["iN"]) inits[paste0("V_",i)] <- init["V"] / (parms["iN"]+1) #add initial values of protected in transit
        }else{
          for(i in 0:parms["iN"]) inits[paste0("V_",i)] <- 0 #add initial values of protected in transit
        }
      }
      inits["R"] <- 0
      if(!is.null(init)) inits["R"] <- init["R"]
      
    }else{
      assign("mod", SPRB_mod)
      inits <- c(S=100, V_0=0)
      if(!is.null(init)) inits[c("S","V_0")] <- init[c("S","V")]
      if(parms["iN"] > 0 & parms["iN"] >= parms["iB"]){
        if(!is.null(init)){
          for(i in 0:parms["iN"]) inits[paste0("V_",i)] <- init["V"] / (parms["iN"]+1) #add initial values of protected in transit
        }else{
          for(i in 0:parms["iN"]) inits[paste0("V_",i)] <- 0 #add initial values of protected in transit
        }
        
      }else if(parms["iN"] > 0 & parms["iN"] < parms["iB"]){
        if(!is.null(init)){
          for(i in 0:parms["iN"]) inits[paste0("V_",i)] <- init["V"] / (parms["iN"]+1) #add initial values of protected in transit
          for(i in (parms["iN"]+1):parms["iB"]) inits[paste0("V_",i)] <- 0 #add initial values of resusceptible in transit
        }else{
          for(i in 0:parms["iB"]) inits[paste0("V_",i)] <- 0 #add initial values of protected/resusceptible in transit
        }
      }
      
      inits["B"] <- 0
      inits["R"] <- 0
      if(!is.null(init)) inits[c("B","R")] <- init[c("B","R")]
    }
    out <- as.data.frame(ode(method=method, inits, times, mod, parms))
    
  }else if(type == "SV"){
    if(parms["m_dur"] == Inf){
      assign("mod", SP_mod) #SP_mod fully corresponds to SV_mod
      inits <- c(S=100, V=0)
      if(!is.null(init)) inits[c("S","V")] <- init[c("S","V")]
      
    }else if(parms["m_dur"] != Inf & parms["boost"] == Inf){
      assign("mod", SVR_mod)
      inits <- c(S=100, V=0, R=0)
      if(!is.null(init)) inits[c("S","V","R")] <- init[c("S","V","R")]
    }else{
      assign("mod", SVR_mod)
      inits <- c(S=100, V=0, R=0)
      if(!is.null(init)) inits[c("S","V","R")] <- init[c("S","V","R")]
    }
    out <- as.data.frame(dede(method=method, inits, times, mod, parms))
  }
}

plotSV <- function(out, dat, i, ...){
    n <- nrow(out[[i]]) #final index (end) of simulated time interval
    if(dat[i,"m_dur"] == Inf){
      plot(out[[i]]$time, out[[i]]$V, ylim=c(0,120), ylab="portion of total cohort (%)", xlab="time (years)", type="l", ...)
      lines(out[[i]]$time, out[[i]]$S, lty=2)

      legend("topright", legend=c("protected","susceptible"), lty=1:2, cex=0.8, horiz=FALSE)

      
      text(x=out[[i]]$time[n]*0.975, y=out[[i]]$V[n]+(0.05*100), labels=format(round(out[[i]]$V[n],1),nsmall=1), cex=0.8)
      text(x=out[[i]]$time[n]*0.975, y=out[[i]]$S[n]-(0.05*100), labels=format(round(out[[i]]$S[n],1),nsmall=1), cex=0.8)
      
      title(main="SV model - lifelong protection")
    }else if(dat[i,"m_dur"] != Inf & dat[i,"boost"] == Inf){
      plot(out[[i]]$time, out[[i]]$V, ylim=c(0,120), ylab="portion of total cohort (%)", xlab="time (years)", type="l", ...)
      lines(out[[i]]$time, out[[i]]$S, lty=2)
      lines(out[[i]]$time, out[[i]]$R, lty=3)

      legend("topright", legend=c("protected","susceptible","resusceptible"), lty=1:3, cex=0.8, horiz=FALSE)
      
      text(x=out[[i]]$time[n]*0.975, y=out[[i]]$V[n]+(0.05*100), labels=format(round(out[[i]]$V[n],1),nsmall=1), cex=0.8)
      text(x=out[[i]]$time[n]*0.975, y=out[[i]]$S[n]-(0.05*100), labels=format(round(out[[i]]$S[n],1),nsmall=1), cex=0.8)
      text(x=out[[i]]$time[n]*0.975, y=out[[i]]$R[n]+(0.05*100), labels=format(round(out[[i]]$R[n],1),nsmall=1), cex=0.8)
      
      title(main=paste("SVR model - protection is lost over time\ntime until protection loss:",dat[i,"m_dur"],"+/-",round(sqrt(dat[i,"var_dur"]),2),"years"))
    }else{
      plot(out[[i]]$time, out[[i]]$V, ylim=c(0,120), ylab="portion of total cohort (%)", xlab="time (years)", type="l", ...)
      lines(out[[i]]$time, out[[i]]$S, lty=2)
      lines(out[[i]]$time, out[[i]]$R, lty=3)

      legend("topright", legend=c("protected","susceptible","resusceptible"), lty=1:3, cex=0.8, horiz=FALSE)
      
      text(x=out[[i]]$time[n]*0.975, y=out[[i]]$V[n]+(0.05*100), labels=format(round(out[[i]]$V[n],1),nsmall=1), cex=0.8)
      text(x=out[[i]]$time[n]*0.975, y=out[[i]]$S[n]-(0.05*100), labels=format(round(out[[i]]$S[n],1),nsmall=1), cex=0.8)
      text(x=out[[i]]$time[n]*0.975, y=out[[i]]$R[n]+(0.05*100), labels=format(round(out[[i]]$R[n],1),nsmall=1), cex=0.8)
      
      title(main=paste("SVRB model - boosting possible after:",dat[i,"boost"],"years\ntime until protection loss:",dat[i,"m_dur"],"+/-",round(sqrt(dat[i,"var_dur"]),2),"years"))
    }
    grid(equilogs = FALSE)
}

plotProt <- function(dat, i, ylab="portion of protected (%) among vaccinated", xlab="time since vaccination (years)", ...){
  if(dat[i,"iN"] != Inf){
    curve(100-100*pgamma(x, shape = dat[i,"m_dur"]*(dat[i,"gamma"]+dat[i,"nu"]), rate = dat[i,"gamma"]+dat[i,"nu"]), ylim=c(0,100), ylab=ylab, xlab=xlab, ...)
    abline(v=dat[i,"m_dur"], lty=2, col="blue")
  }else{
    plot(c(1,40),c(0,100), type="n", ylab=ylab, xlab=xlab, ...)
    abline(h=100)
  }
  grid(equilogs = FALSE)
}
plotBoost <- function(dat, i, ylab="portion of boosted (%) among vaccinated", xlab="time since vaccination (years)", ...){
  if(dat[i,"iN"] != Inf){
    curve(100*pexp(x-dat[i,"boost"], rate = dat[i,"delta"]), ylim=c(0,100), ylab=ylab, xlab=xlab, ...)
    abline(v=dat[i,"boost"], lty=2, col="blue")
  }else{
    plot(c(1,40),c(0,100), type="n", ylab=ylab, xlab=xlab, ...)
    abline(h=0)
  }
  grid(equilogs = FALSE)
}

PortionRate <- Vectorize(function(portion, years=10){
  if(portion <= 0 | portion >= 1 | years == Inf | years <= 0){
    tmp <- 0
  }else{
    tmp <- optimize(function(x,portion,years) (qexp(portion,x)-years)^2, lower = 0, upper = 10, portion=portion, years=years)$minimum
  }
  return(tmp)
}, vectorize.args = c("portion","years"))

convert_input2dat <- function(input, cached=F){
  if(cached){
    UserPars <- c("deathPortion"=input$deathPortion_C/100,
                  "vaccPortion"=input$vaccPortion_C/100,
                  "m_dur"=as.numeric(input$m_dur_C),
                  "cv_dur"=input$cv_dur_C/100,
                  "boost"=as.numeric(input$boost_C),
                  "boostPortion"=input$boostPortion_C/100)
    UserPars[6] <- ifelse(UserPars[3] != Inf & UserPars[5] == Inf, unique(dat$boostPortion)[1], UserPars[6])
    UserPars[5] <- ifelse(UserPars[3] == Inf, unique(dat$boost)[1], UserPars[5])
    UserPars[6] <- ifelse(UserPars[3] == Inf, unique(dat$boostPortion)[1], UserPars[6])
    
    which(apply(dat[,1:6],1, function(x) sum(x == UserPars) == 6))
  }else{
    UserPars <- c("deathPortion"=input$deathPortion_S/100,
                  "vaccPortion"=input$vaccPortion_S/100,
                  "m_dur"=as.numeric(input$m_dur_S),
                  "cv_dur"=input$cv_dur_S/100,
                  "boost"=as.numeric(input$boost_S),
                  "boostPortion"=input$boostPortion_S/100,
                  "f_nu"=input$vaccdeathRR_S)
    
    tab <- NULL
    tab["mu"] <- tab["nu"] <- PortionRate(UserPars["deathPortion"], 1) #determine and constrain birth rate to mortality rate
    tab["var_dur"] <- (UserPars["cv_dur"]*UserPars["m_dur"])^2 #variance of protection duration (years)
    tab["lambda"] <- PortionRate(UserPars["vaccPortion"], 1) #vaccination rate per year
    tab["gamma"] <- ifelse(UserPars["m_dur"] != Inf, UserPars["m_dur"]/tab["var_dur"] - tab["nu"], 0) #transition rate per year
    tab["delta"] <- PortionRate(UserPars["boostPortion"], 1) #boost rate per year
    tab["iN"] <- round(tab["var_dur"]*(tab["gamma"]+tab["nu"])^2 -1) #total number of transition compartments
    tab["iB"] <- round((tab["var_dur"]*(tab["gamma"]+tab["nu"])^2 -1)*(UserPars["boost"]/UserPars["m_dur"])) #first boost transition compartment
    tab["iB"][is.na(tab["iB"])] <- tab["iN"][is.na(tab["iB"])]
    
    data.frame(matrix(c(UserPars,tab), nrow=1, dimnames = list(NULL,c(names(UserPars),names(tab))) ))  
  }
}

seqb <- function(x, length = 50, digits = 1){
  if(any(x == Inf)){
    x <- x[x != Inf]
    rx <- log(range(x))
    ret <- c(seq(from = rx[1], to = rx[2], length = length-1),Inf)
  }else{
    rx <- log(range(x))
    ret <- seq(from = rx[1], to = rx[2], length = length)
  }
  return(round(exp(ret),digits))
}