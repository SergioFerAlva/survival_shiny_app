library(shiny)
library(survival)
library(survminer)
library(dplyr)

# Load dataset safely
lung <- survival::lung

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
      h3("Descriptive Survival Statistics"),
      tableOutput("surv_stats"),
      
      h3("Kaplan-Meier Curve"),
      plotOutput("km_plot")
    )
  )
)

server <- function(input, output) {
  
  # Reactive survival model
  surv_object <- reactive({
    if (input$strata == "None") {
      survfit(Surv(time, status) ~ 1, data = lung)
    } else {
      survfit(Surv(time, status) ~ sex, data = lung)
    }
  })
  
  # Kaplan-Meier Plot
  output$km_plot <- renderPlot({
    
    fit <- surv_object()
    
    p <- ggsurvplot(
      fit,
      data = lung,
      risk.table = FALSE,
      pval = input$strata != "None",
      conf.int = TRUE
    )
    
    print(p$plot)
  })
  
  # Descriptive Survival Statistics
  output$surv_stats <- renderTable({
    
    fit <- surv_object()
    
    # Get summary table
    fit_table <- summary(fit)$table
    
    # ---- CASE 1: No stratification (vector) ----
    if (is.null(dim(fit_table))) {
      
      total_n <- nrow(lung)
      total_events <- sum(lung$status)
      
      med <- as.numeric(fit_table["median"])
      
      surv_365 <- summary(fit, times = 365, extend = TRUE)$surv[1]
      
      data.frame(
        Group = "Overall",
        Total_Patients = total_n,
        Events = total_events,
        Event_Rate_Percent = round(100 * total_events / total_n, 1),
        Median_Survival_Time = round(med, 1),
        Survival_Probability_365_Days = round(surv_365, 3)
      )
      
    } else {
      
      # ---- CASE 2: Stratified (matrix) ----
      
      summary_df <- lung %>%
        group_by(sex) %>%
        summarise(
          Total_Patients = n(),
          Events = sum(status),
          .groups = "drop"
        )
      
      med <- as.numeric(fit_table[, "median"])
      surv_365 <- summary(fit, times = 365, extend = TRUE)$surv
      
      summary_df$Event_Rate_Percent <- round(
        100 * summary_df$Events / summary_df$Total_Patients, 1
      )
      
      summary_df$Median_Survival_Time <- round(med, 1)
      summary_df$Survival_Probability_365_Days <- round(surv_365, 3)
      
      summary_df
    }
    
  })
  }

shinyApp(ui = ui, server = server)