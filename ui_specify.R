fixedRow(
  column(width = 4, #sidebar panel
         h5(em("Implied model parameters")),
         wellPanel(
           tableOutput("ImplParmsS")
         ),
         wellPanel(
           sliderInput("rateconv", label = "Crude annual rate", post = "%",
                       min=1, max=99, step = 1, value = 1),
           # sliderTextInput("ratetime", label = "Duration", post = " years", grid = TRUE,
           #              choices = c(0.5, 1, 2, 5, 10, 20), selected = 1),
           plotOutput("grateconv")
         )
  ),
  column(width = 8, #main panel
         fixedRow( #plot output subpanel
           column(width = 12, 
                  h5(em("R model expression to be evaluated by eval()")),
                  verbatimTextOutput("syntax")
           )
         )
  )
)