fixedRow(useShinyjs(),
  column(width = 4, #sidebar panel
         wellPanel(
           sliderTextInput("deathPortion_S",
                           "Crude annual mortality rate",
                           post = " %",
                           choices = seqb(c(0.25,32), digits = 2),
                           selected = 1.35),
           sliderTextInput("vaccPortion_S",
                           "Crude annual vaccination rate",
                           post = " %",
                           choices = c(0.25, seq(0.5,5, by=0.5), seq(1,96, by=1)),
                           selected=20),
           sliderTextInput("vaccEffect_S",
                           "Vaccination efficacy",
                           post = "",
                           choices = seq(0.25,1, 0.05),
                           selected=1),
           sliderTextInput("m_dur_S",
                           "Mean duration until protection loss",
                           post = " years",
                           choices = seqb(unique(dat$m_dur))),
           sliderTextInput("cv_dur_S",
                           "CV of duration until protection loss",
                           post = " %",
                           choices = seqb(unique(dat$cv_dur)*100)),
           sliderTextInput("boost_S",
                           "Minimal boosting interval",
                           post = " years",
                           choices = seqb(c(1,40,Inf)),
                           selected=8),
           sliderTextInput("boostPortion_S",
                           "Crude annual boosting rate",
                           post = " %",
                           choices = c(0.25, seq(0.5,5, by=0.5), seq(1,96, by=1)),
                           selected=30),
           br(),
           h5(em("Cohort state at 0 years")),
           sliderTextInput("V_initial",
                           "Portion of protected among total",
                           post = " %",
                           choices = seq(0,100,2),
                           selected = 0),
           sliderTextInput("SR_initial",
                           "Portion of resusceptible among unprotected",
                           post = " %",
                           choices = seq(0,100,2),
                           selected = 0)
         )
  ),
  column(width = 8, #main panel
        fixedRow( #plot output subpanel
          column(width = 12, 
                 tabsetPanel(
                   tabPanel("Coverage", plotOutput("VaccPlotS0", click = "locator")),
                   tabPanel("Waning", plotOutput("VaccPlotS1", click = "locator")),
                   tabPanel("Boosting", plotOutput("VaccPlotS2", click = "locator"))
                 ), br(),br(),br(),br(),br(),br()
          )
        ),
        fixedRow( #plot control subpanel
          column(width = 5, wellPanel(tableOutput("CursorPos")), align="center"),
          column(width = 7, align ="center",
                wellPanel(
                  downloadButton('downloadData', 'Download Data'),
                  downloadButton('downloadPlot', 'Download Plot') 
                )
          )
        )
  )
)