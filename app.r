library(shiny)
library(survival)
library(survminer)
library(dplyr)

# Load dataset
data(lung)

# Clean data
lung <- lung %>%
  mutate(
    status = as.numeric(status == 2),
    sex = factor(sex, labels = c("Male", "Female"))
  )

ui <- fluidPage(
  titlePanel("Interactive Survival Explorer"),

  sidebarLayout(
    sidebarPanel(
      selectInput("strata",
                  "Stratify by:",
                  choices = c("None", "sex"),
                  selected = "None")
    ),

    mainPanel(
      h3("Dataset Summary"),
      verbatimTextOutput("summary"),

      h3("Kaplan-Meier Curve"),
      plotOutput("km_plot")
    )
  )
)

server <- function(input, output) {

  output$summary <- renderPrint({
    summary(lung)
  })

  output$km_plot <- renderPlot({

    if (input$strata == "None") {
      fit <- survfit(Surv(time, status) ~ 1, data = lung)
    } else {
      formula <- as.formula(paste("Surv(time, status) ~", input$strata))
      fit <- survfit(formula, data = lung)
    }

    survminer::ggsurvplot(
      fit,
      data = lung,
      risk.table = FALSE,  # temporarily disable
      pval = TRUE,
      conf.int = TRUE
    )
  })
}

shinyApp(ui = ui, server = server)