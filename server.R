library(shiny)
library(shinyWidgets)

source("SVpipe.R")
source("SV.R")
source("SVR.R")
source("SVRB.R")

#dat <- readRDS("SVdat.rds")
#out <- readRDS("SVsim.rds")

function(input, output) {

output$ImplParmsC <- renderTable({
    i <- convert_input2dat(input, cached = T)
    
    if(dat[i,"m_dur"] == Inf){
      tab <- dat[i,c("nu","lambda")]
      data.frame("Parameter"= c(HTML("&nu;"),HTML("&lambda;")),
                 "Value"=as.numeric(tab))
    }else if(dat[i,"m_dur"] != Inf & dat[i,"boost"] == Inf){
      tab <- dat[i,c("nu","lambda","gamma","iN")]
      data.frame("Parameter"= c(HTML("&nu;"),HTML("&lambda;"),HTML("&gamma;"),"iN"),
                 "Value"=as.numeric(tab))
    }else if(dat[i,"iB"] > dat[i,"iN"]){
      tab <- dat[i,c("nu","lambda","gamma","iN","delta","iB")]
      data.frame("Parameter"= c(HTML("&nu;"),HTML("&lambda;"),HTML("&gamma;"),"iN",HTML("&delta;"),"iB"),
                 "Value"=as.numeric(tab))
    }else{
      tab <- dat[i,c("nu","lambda","gamma","iN","delta","iB")]
      data.frame("Parameter"= c(HTML("&nu;"),HTML("&lambda;"),HTML("&gamma;"),"iN",HTML("&delta;"),"iB"),
                 "Value"=as.numeric(tab))
    }
  }, digits = 3, sanitize.text.function = function(x) x)
  
  #generate plot from cached simulations
  output$VaccPlotC <- renderPlot({
    i <- convert_input2dat(input, cached = T)
    
    layout(matrix(c(1,2,1,3), nrow=2, ncol=2))
    plotSV(out, dat, i, xlim=c(0,100))
    plotProt(dat, i, xlim=c(1,60), log="x")
    plotBoost(dat, i, xlim=c(1,60), log="x")
  }, width = 700, height = 800)

  
  #perform simulation
  simCoverage <- reactive({
    dat <- convert_input2dat(input, cached = F)

    times <- seq(0, 100, length=201)
    parms <- with(dat[1,], c(mu=mu, nu=nu, lambda=lambda, delta=delta, gamma=gamma, iN=iN, iB=iB, m_dur=m_dur, boost=boost, f_nu=f_nu))
    out <- list(SVpipe(times, parms, init = c("S"=(100-input$V_initial)*(1-input$SR_initial/100),
                                              "V"=input$V_initial,
                                              "R"=(100-input$V_initial)*(input$SR_initial/100),
                                              "B"=0)
                       )
                )
    out <- PruneOutput(dat,out)
    return(list("out"=out, "dat"=dat))
  })
  
  
  #generate plots
  output$VaccPlotS0 <- renderPlot({
    res <- simCoverage()
    plotSV(res[["out"]], res[["dat"]], 1, xlim=c(0,100))
  }, width = "auto", height = 500)
  
  output$VaccPlotS1 <- renderPlot({
    dat <- convert_input2dat(input, cached = F)
    plotProt(dat, 1, xlim=c(1,60), log="x")
  }, width = "auto", height = 500)
  
  output$VaccPlotS2 <- renderPlot({
    dat <- convert_input2dat(input, cached = F)
    plotBoost(dat, 1, xlim=c(1,60), log="x")
  }, width = "auto", height = 500)
  
  #ui input monitoring
  observe({
    toggleState("cv_dur_C", condition = input$m_dur_C != Inf)
    toggleState("boost_C", condition = input$m_dur_C != Inf)
    toggleState("boostPortion_C", condition = input$m_dur_C != Inf & input$boost_C != Inf)
    
    toggleState("cv_dur_S", condition = input$m_dur_S != Inf)
    toggleState("boost_S", condition = input$m_dur_S != Inf)
    toggleState("SR_initial", condition = input$m_dur_S != Inf)
    toggleState("boostPortion_S", condition = input$m_dur_S != Inf & input$boost_S != Inf)
  })
  
  #ui utility functions
  output$CursorPos <- renderTable({
    if(is.null(input$locator)){
      data.frame("Time"=0, "Percent"=0)
    }else{
      data.frame("Time"=input$locator$x, "Percent"=input$locator$y)
    }
  }, digits = 1)
  
  output$downloadData <- downloadHandler(
    filename = function() { paste('coverage.csv', sep='') },
    content = function(file) {
      res <- simCoverage()
      write.csv(res[["out"]], file, row.names = FALSE)
    }

  )
  output$downloadPlot <- downloadHandler(
    filename = function() { paste('coverage.pdf', sep='') },
    content = function(file) {
      pdf(file)
      res <- simCoverage()
      plotSV(res[["out"]], res[["dat"]], 1, xlim=c(0,100))
      dev.off()
    }
  )
  
  #generate table of implied model parameters
  output$ImplParmsS <- renderTable({
    dat <- convert_input2dat(input, cached = F)
    i <- 1
    
    if(dat[i,"m_dur"] == Inf){
      tab <- dat[i,c("nu","lambda")]
      data.frame("Parameter"= c(HTML("&nu;"),HTML("&lambda;")),
                 "Value"=as.numeric(tab))
    }else if(dat[i,"m_dur"] != Inf & dat[i,"boost"] == Inf){
      tab <- dat[i,c("nu","lambda","gamma","iN")]
      data.frame("Parameter"= c(HTML("&nu;"),HTML("&lambda;"),HTML("&gamma;"),"iN"),
                 "Value"=as.numeric(tab))
    }else if(dat[i,"iB"] > dat[i,"iN"]){
      tab <- dat[i,c("nu","lambda","gamma","iN","delta","iB")]
      data.frame("Parameter"= c(HTML("&nu;"),HTML("&lambda;"),HTML("&gamma;"),"iN",HTML("&delta;"),"iB"),
                 "Value"=as.numeric(tab))
    }else{
      tab <- dat[i,c("nu","lambda","gamma","iN","delta","iB")]
      data.frame("Parameter"= c(HTML("&nu;"),HTML("&lambda;"),HTML("&gamma;"),"iN",HTML("&delta;"),"iB"),
                 "Value"=as.numeric(tab))
    }
  }, digits = 3, sanitize.text.function = function(x) x)
  
  #plot of rate conversion
  output$grateconv <- renderPlot({
    p <- seqb(c(0.01,0.99), length = 50, digits = 3)
    
    plot(p,PortionRate(p,1), type="l", xlab="crude annual rate", ylab="instantaneous annual rate", log="xy")
    grid(equilogs = F)
    lines(p,p, lty=2, col="darkgrey")
    
    coords <- round(c(input$rateconv/100, PortionRate(input$rateconv/100,1)),3)
    points(coords[1], coords[2], pch=19)
    text(0.05,2, labels = paste("rate =", coords[2]))
  }, bg="transparent")
  
  
  #generate output of model syntax
  output$syntax <- renderPrint({
    dat <- convert_input2dat(input, cached = F)
    i <- 1
    
    if(dat[i,"m_dur"] == Inf){
      mod <- readLines("SV.R", n = 24)
      env1 <- c('times <- seq(0, 100, length=101)',
               with(dat[i,], paste('parms <- c("nu"=',nu,',"lambda"=',lambda,',"f_nu"=',f_nu,')')),
               paste0('inits <- c("S"=',(100-input$V_initial)*(1-input$SR_initial/100),', "V"=',input$V_initial,')'),
               'assign("mod", SV_mod)')
      env3 <- c('V <- out$V', 'R <- rep(0, length(out$time))')
    }else if(dat[i,"m_dur"] != Inf & dat[i,"boost"] == Inf){
      mod <- readLines("SVR.R", n = 41)
      env1 <- c('times <- seq(0, 100, length=101)',
               with(dat[i,], paste('parms <- c("nu"=',nu,',"lambda"=',lambda,',"gamma"=',gamma,',"iN"=',iN,',"f_nu"=',f_nu,')')),
               paste0('init <- c("S"=',(100-input$V_initial)*(1-input$SR_initial/100),', "V"=',input$V_initial,', "R"=',(100-input$V_initial)*(input$SR_initial/100),')'),
               'assign("mod", SVR_mod)',
               'inits <- NULL',
               'inits[c("S","V_0")] <- init[c("S","V")]',
               'if(parms["iN"] > 0){',
               '  for(i in 1:parms["iN"]) inits[paste0("V_",i)] <- 0 #add initial values of protected in transit',
               '}',
               'inits["R"] <- init["R"]'
              )
      env3 <- c('V <- apply(out[3:(ncol(out)-1)],1,sum)', 'R <- out$R')
    }else{
      mod <- readLines("SVRB.R", n = 72)
      env1 <- c('times <- seq(0, 100, length=101)',
               with(dat[i,],paste('parms <- c("nu"=',nu,',"lambda"=',lambda,',"delta"=',delta,',"gamma"=',gamma,',"iN"=',iN,',"iB"=',iB,',"f_nu"=',f_nu,')')),
               paste0('init <- c("S"=',(100-input$V_initial)*(1-input$SR_initial/100),', "V"=',input$V_initial,', "R"=',(100-input$V_initial)*(input$SR_initial/100),', "B"=0)'),
               'assign("mod", SVRB_mod)',
               'inits <- NULL',
               'inits <- c(S=init["S"], V_0=init["V"])',
               'if(parms["iN"] > 0 & parms["iN"] >= parms["iB"]){',
               '  for(i in 1:parms["iN"]) inits[paste0("V_",i)] <- 0 #add initial values of protected in transit',
               '}else if(parms["iN"] > 0 & parms["iN"] < parms["iB"]){',
               '  for(i in 1:parms["iB"]) inits[paste0("V_",i)] <- 0 #add initial values of protected/resuceptible in transit',
               '}',
               'inits[c("B","R")] <- init[c("B","R")]'
               )
      if(dat[i,"iN"] >= dat[i,"iB"]){
        env3 <- c('V <- apply(out[3:(ncol(out)-2)],1,sum)', 'R <- out$R')
      }else if(dat[i,"iN"] < dat[i,"iB"]){
        env3 <- c('V <- apply(out[3:(ncol(out)-2-parms["iB"]+parms["iN"])],1,sum)', 'R <- apply(out[(ncol(out)-1-parms["iB"]+parms["iN"]):(ncol(out)-2)],1,sum) + out$R')
      }
    }
    env2 <- c('out <- as.data.frame(ode(method="lsoda", inits, times, mod, parms))', 'S <- out$S')
    env4 <- c('out <- data.frame("time"=out$time, "S"=S, "V"=V, "R"=R)',
              'suppressWarnings(rm(times,init,inits,S,V,R,i,mod))',
              'matplot(out[,1],out[,-1], type="l", ylab="portion of total cohort (%)", xlab="time (years)", lty=c(3,1,2))')
    
    temp <- c(mod,env1,env2,env3,env4)
    print(parse(text = temp))
  })
}
