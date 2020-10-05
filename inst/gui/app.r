#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(shinythemes)
library(DT)
# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  withMathJax(),
  # Application title
  titlePanel("Contional Prevalence Estimation"),
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      h4("Entre Data:"),
      numericInput("R1",
                   withMathJax('\\(R_1 = \\)'),
                   value = 30),
      numericInput("R2",
                   withMathJax('\\(R_2 = \\)'),
                   value = 0),
      numericInput("R3",
                   withMathJax('\\(R_3 = \\)'),
                   value = 10),
      numericInput("R4",
                   withMathJax('\\(R_4 = \\)'),
                   value = 960),
      numericInput("pi0",
                   withMathJax('\\(\\pi_0 = \\)'),
                   value = 0.01),
      HTML("<br>"),
      h4("Summary Tab"),
      selectInput("select", "Select Estimator",
                  choices = list("Condition MLE" = 1, "Marginal MLE" = 2,
                                 "Method of Moments" = 3, "Survey MLE" = 4), selected = 1),
      numericInput("gamma",
                   withMathJax('\\(\\gamma = \\)'),
                   value = 0.05),
      HTML("<br>"),
      h4("Misclassification error:"),
      numericInput("alpha0",
                   withMathJax('\\(\\alpha_0 = \\)'),
                   value = 0),
      numericInput("alpha",
                   withMathJax('\\(\\alpha = \\)'),
                   value = 0),
      numericInput("beta0",
                   withMathJax('\\(\\beta_0 = \\)'),
                   value = 0),
      numericInput("beta",
                   withMathJax('\\(\\beta = \\)'),
                   value = 0)
    ),
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Comparison", plotOutput("plot", width = "100%")),
                  tabPanel("Summary", verbatimTextOutput("summary")),
                  tabPanel("Table", DT::dataTableOutput("table"))
      ))
  )
)
# Define server logic required to draw a histogram
server <- function(input, output) {
  output$plot <- renderPlot({
    cmle80 = conditional_mle(R1 = input$R1, R2 = input$R2, R3 = input$R3, R4 = input$R4, pi0 = input$pi0, gamma = 0.2)
    cmle95 = conditional_mle(R1 = input$R1, R2 = input$R2, R3 = input$R3, R4 = input$R4, pi0 = input$pi0, gamma = 0.05)
    cmle99 = conditional_mle(R1 = input$R1, R2 = input$R2, R3 = input$R3, R4 = input$R4, pi0 = input$pi0, gamma = 0.01)
    mmle80 = marginal_mle(R1 = input$R1, R3 = input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, pi0 = input$pi0, gamma = 0.2)
    mmle95 = marginal_mle(R1 = input$R1, R3 = input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, pi0 = input$pi0, gamma = 0.05)
    mmle99 = marginal_mle(R1 = input$R1, R3 = input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, pi0 = input$pi0, gamma = 0.01)
    mom80 = moment_estimator(R3 = input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, pi0 = input$pi0, gamma = 0.2)
    mom95 = moment_estimator(R3 = input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, pi0 = input$pi0, gamma = 0.05)
    mom99 = moment_estimator(R3 = input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, pi0 = input$pi0, gamma = 0.01)
    sur80 = survey_mle(R = input$R1 + input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, gamma = 0.2)
    sur95 = survey_mle(R = input$R1 + input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, gamma = 0.05)
    sur99 = survey_mle(R = input$R1 + input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, gamma = 0.01)
    cols = c("#adc2eb", "#7094db", "#2e5cb8", "#e68a00")
    cmle = c(cmle80$estimate, cmle80$ci_asym, cmle95$ci_asym, cmle99$ci_asym)
    mmle = c(mmle80$estimate, mmle80$ci_asym, mmle95$ci_asym, mmle99$ci_asym)
    mom_asy = c(mom80$estimate, mom80$ci_asym, mom95$ci_asym, mom99$ci_asym)
    mom_cp = c(mom80$estimate, mom80$ci_cp, mom95$ci_cp, mom99$ci_cp)
    sur_asy = c(sur80$estimate, sur80$ci_asym, sur95$ci_asym, sur99$ci_asym)
    sur_cp = c(sur80$estimate, sur80$ci_cp, sur95$ci_cp, sur99$ci_cp)
    mat = cbind(cmle, mmle, mom_asy, mom_cp, sur_asy, sur_cp)
    par(lend=1, mar = c(5,8,2,4))
    plot(NA, xlim = c(min(mat) - 0.1*diff(range(mat)), max(mat) + 0.1*diff(range(mat))), ylim = c(0.5, 6.5), axes = FALSE, ann = FALSE)
    grid()
    box()
    axis(1, cex = 1.2)
    mtext("Estimated Prevalence", side = 1, line = 3, cex = 1.2)
    delta = 0.06
    meth = c("Cond. MLE (A)", "Marg. MLE (A)", "MoM (A)", "MoM (CP)", "Sur. MLE (A)", "Sur. MLE (CP)")
    axis(2, at = 1:6, labels = meth, las = 2, cex = 1.2)
    for (i in 1:6){
      lines(c(mat[6,i], mat[7,i]), c(i, i), lwd = 2, col = cols[1])
      lines(c(mat[6,i], mat[6,i]), c(i - delta, i + delta), lwd = 2, col = cols[1])
      lines(c(mat[7,i], mat[7,i]), c(i - delta, i + delta), lwd = 2, col = cols[1])
      lines(c(mat[4,i], mat[5,i]), c(i, i), lwd = 5, col = cols[2])
      lines(c(mat[4,i], mat[4,i]), c(i - delta, i + delta), lwd = 2, col = cols[2])
      lines(c(mat[5,i], mat[5,i]), c(i - delta, i + delta), lwd = 2, col = cols[2])
      lines(c(mat[2,i], mat[3,i]), c(i, i), lwd = 8, col = cols[3])
      lines(c(mat[2,i], mat[2,i]), c(i - delta, i + delta), lwd = 2, col = cols[3])
      lines(c(mat[3,i], mat[3,i]), c(i - delta, i + delta), lwd = 2, col = cols[3])
      points(mat[1,i], i, col = cols[4], pch =16, cex = 2)
    }
    legend("bottomright", c("Estimate", "Confidence level:", "80%", "95%", "99%"), bty = "n",
           pch = c(16, NA, NA, NA, NA), lwd = c(NA, NA, 8, 5, 2),
           col = c(cols[4], NA, cols[3:1]), pt.cex = c(2, NA, NA, NA, NA), cex = 1.05)
  }, height = 500, width = 600)
  output$summary <- renderPrint({
    if (input$select == 1){
      conditional_mle(R1 = input$R1, R2 = input$R2, R3 = input$R3, R4 = input$R4, pi0 = input$pi0, gamma = input$gamma)
    }else{
      if (input$select == 2){
        marginal_mle(R1 = input$R1, R3 = input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, pi0 = input$pi0, gamma = input$gamma)
      }else{
        if (input$select == 3){
          moment_estimator(R3 = input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, pi0 = input$pi0, gamma = input$gamma)
        }else{
          survey_mle(R = input$R1 + input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, gamma = input$gamma)
        }
      }
    }
  })
  output$table <- DT::renderDataTable(DT::datatable({
    cmle = conditional_mle(R1 = input$R1, R2 = input$R2, R3 = input$R3, R4 = input$R4, pi0 = input$pi0, gamma = input$gamma)
    mmle = marginal_mle(R1 = input$R1, R3 = input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, pi0 = input$pi0, gamma = input$gamma)
    mom = moment_estimator(R3 = input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, pi0 = input$pi0, gamma = input$gamma)
    sur = survey_mle(R = input$R1 + input$R3, n = input$R1 + input$R2 + input$R3 + input$R4, gamma = input$gamma)
    cmle = c(cmle$estimate, cmle$ci_asym)
    mmle = c(mmle$estimate, mmle$ci_asym)
    mom_asy = c(mom$estimate, mom$ci_asym)
    mom_cp = c(mom$estimate, mom$ci_cp)
    sur_asy = c(sur$estimate, sur$ci_asym)
    sur_cp = c(sur$estimate, sur$ci_cp)
    A = round(rbind(cmle, mmle, mom_asy, mom_cp, sur_asy, sur_cp), 6)
    colnames(A) = c("Estimate", "CI (lower bound)", "CI (upper bound)")
    rownames(A) = c("Conditional MLE (Asy.)", "Marginal MLE (Asy.)", "Method of Moments (Asy.)",
                    "Method of Moments (CP)", "Survey MLE (Aym.)", "Survey MLE (CP)")
    A
  }))
}
# Run the application
shinyApp(ui = ui, server = server)
