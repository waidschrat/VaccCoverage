library(shiny)
library(shinyWidgets)
library(shinyjs)

dat <- readRDS("SVdat.rds")
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

navbarPage("Vaccination Coverage", selected = "Simulation", id="mode", fluid = FALSE, collapsible = TRUE,
           navbarMenu("Mode", icon = icon("code"),
                      tabPanel("Simulation", source("ui_simulate.R", local = TRUE)$value),
                      #tabPanel("Cached", source("ui_cached.R", local = TRUE)),
                      tabPanel("Specification", source("ui_specify.R", local = TRUE)$value)
                      )
)