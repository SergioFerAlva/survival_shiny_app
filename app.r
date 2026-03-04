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
      plotOutput("km_plot"),
      
      h3("Cox Proportional Hazards Model"),
      tableOutput("cox_results"),
      verbatimTextOutput("cox_interpretation"),
      
      h4("Multivariable Cox Model (Adjusted)"),
      tableOutput("cox_multi_results"),
      plotOutput("forest_plot"),
      
      h3("Proportional Hazards Assumption Test"),
      tableOutput("ph_test"),
      verbatimTextOutput("ph_interpretation"),
      
      h4("Schoenfeld Residual Diagnostics"),
      plotOutput("schoenfeld_plot")
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
      conf.int = TRUE,
      xlab = "Time (days)",
      ylab = "Survival Probability",
      ggtheme = theme_minimal()
      
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
  
  # Cox proportional hazards model
  output$cox_results <- renderTable({
    
    if (input$strata != "sex") {
      return(data.frame(
        Message = "Cox model available when stratifying by sex."
      ))
    }
    
    model <- coxph(Surv(time, status) ~ sex, data = lung)
    model_summary <- summary(model)
    
    hr <- model_summary$coefficients[,"exp(coef)"]
    lower_ci <- model_summary$conf.int[,"lower .95"]
    upper_ci <- model_summary$conf.int[,"upper .95"]
    p_value <- model_summary$coefficients[,"Pr(>|z|)"]
    
    data.frame(
      Variable = "Female vs Male",
      Hazard_Ratio = round(hr, 3),
      CI_95 = paste0(round(lower_ci,3), " - ", round(upper_ci,3)),
      P_Value = signif(p_value, 3)
    )
    
  })
  
  # Cox interpretation blurb
  output$cox_interpretation <- renderPrint({
    
    if (input$strata != "sex") {
      return(NULL)
    }
    
    model <- coxph(Surv(time, status) ~ sex, data = lung)
    model_summary <- summary(model)
    
    hr <- model_summary$coefficients[,"exp(coef)"]
    p_value <- model_summary$coefficients[,"Pr(>|z|)"]
    
    if (p_value < 0.05) {
      significance <- "statistically significant"
    } else {
      significance <- "not statistically significant"
    }
    
    cat(
      "Interpretation:\n\n",
      "Female patients have a hazard ratio of", round(hr,3),
      "relative to male patients.\n",
      "This difference is", significance,
      "(p =", signif(p_value,3), ")."
    )
  })
  
  # Multivariate Cox
  output$cox_multi_results <- renderTable({
    
    if (input$strata != "sex") {
      return(NULL)
    }
    
    model <- coxph(Surv(time, status) ~ sex + age + ph.ecog, data = lung)
    summary_model <- summary(model)
    
    hr <- summary_model$coefficients[, "exp(coef)"]
    lower <- summary_model$conf.int[, "lower .95"]
    upper <- summary_model$conf.int[, "upper .95"]
    pval <- summary_model$coefficients[, "Pr(>|z|)"]
    
    data.frame(
      Variable = rownames(summary_model$coefficients),
      Hazard_Ratio = round(hr, 3),
      CI_95 = paste0(round(lower,3), " - ", round(upper,3)),
      P_Value = signif(pval, 3)
    )
  })
  
  # Forest plot
    output$forest_plot <- renderPlot({
    
    if (input$strata != "sex") {
      return(NULL)
    }
    
    model <- coxph(Surv(time, status) ~ sex + age + ph.ecog, data = lung)
    
    ggforest(model, data = lung)
  })
  
  # PH Test Table
  output$ph_test <- renderTable({
    
    if (input$strata != "sex") {
      return(NULL)
    }
    
    model <- coxph(Surv(time, status) ~ sex + age + ph.ecog, data = lung)
    
    ph_test <- cox.zph(model)
    
    data.frame(
      Variable = rownames(ph_test$table),
      Chi_Square = round(ph_test$table[, "chisq"], 3),
      P_Value = signif(ph_test$table[, "p"], 3)
    )
  })
  
  # Schoenfeld residual plots
  output$schoenfeld_plot <- renderPlot({
    
    if (input$strata != "sex") {
      return(NULL)
    }
    
    model <- coxph(Surv(time, status) ~ sex + age + ph.ecog, data = lung)
    ph_test <- cox.zph(model)
    
    plot(ph_test)
  })
  
  # PH interpretation
  output$ph_interpretation <- renderPrint({
    
    if (input$strata != "sex") {
      return(NULL)
    }
    
    model <- coxph(Surv(time, status) ~ sex + age + ph.ecog, data = lung)
    ph_test <- cox.zph(model)
    
    pvals <- ph_test$table[, "p"]
    global_p <- pvals[length(pvals)]
    
    if (global_p < 0.05) {
      cat(
        "Interpretation:\n\n",
        "The global Schoenfeld test indicates a violation of the proportional hazards assumption",
        "(p =", signif(global_p,3), ").\n",
        "This suggests that at least one covariate may have time-dependent effects.",
        "Inspection of the Schoenfeld residual plots is recommended to identify which variable violates the assumption."
      )
    } else {
      cat(
        "Interpretation:\n\n",
        "The global Schoenfeld test does not indicate violation of the proportional hazards assumption",
        "(p =", signif(global_p,3), ").\n",
        "The proportional hazards assumption appears to be satisfied for this Cox model."
      )
    }
  })
  
  }

shinyApp(ui = ui, server = server)