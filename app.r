###############################################################
# Interactive Survival Explorer
# Bioinformatics / Biostatistics Shiny portfolio project
###############################################################

library(shiny)
library(survival)
library(survminer)
library(dplyr)

###############################################################
# DEFAULT DATASET
###############################################################

# Load lung dataset from survival package
lung <- survival::lung

# Clean variables
lung <- lung %>%
  mutate(
    status = as.numeric(status == 2),        # convert event indicator to 1/0
    sex = factor(sex, labels = c("Male","Female"))
  )

###############################################################
# USER INTERFACE
###############################################################

ui <- fluidPage(
  
  titlePanel("Interactive Survival Explorer"),
  
  sidebarLayout(
    
    ###########################################################
    # SIDEBAR
    ###########################################################
    
    sidebarPanel(
      
      # Optional CSV upload
      fileInput(
        "file",
        "Upload CSV dataset",
        accept = ".csv"
      ),
      
      # Dynamic UI for variable selection
      uiOutput("var_select"),
      
      # Covariate selection for multivariable Cox model
      selectizeInput(
        "covariates",
        "Covariates for multivariable Cox model",
        choices = NULL,
        multiple = TRUE
      )
      
    ),
    
    ###########################################################
    # MAIN PANEL
    ###########################################################
    
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
      
      h4("Proportional Hazards Diagnostics"),
      tableOutput("ph_test"),
      plotOutput("ph_plot"),
      verbatimTextOutput("ph_interpretation")
      
    )
  )
)

###############################################################
# SERVER LOGIC
###############################################################

server <- function(input, output, session) {
  
  #############################################################
  # DATA INPUT
  #############################################################
  
  # Reactive dataset (default = lung)
  data_input <- reactive({
    
    if (is.null(input$file)) {
      return(lung)
    }
    
    read.csv(input$file$datapath)
    
  })
  
  
  #############################################################
  # VARIABLE SELECTION UI
  #############################################################
  
  output$var_select <- renderUI({
    
    df <- data_input()
    
    tagList(
      
      selectInput(
        "time_var",
        "Time variable:",
        choices = names(df),
        selected = if("time" %in% names(df)) "time" else names(df)[1]
      ),
      
      selectInput(
        "event_var",
        "Event variable (1=event, 0=censored):",
        choices = names(df),
        selected = if("status" %in% names(df)) "status" else names(df)[2]
      ),
      
      selectInput(
        "group_var",
        "Grouping variable:",
        choices = c("None", names(df)),
        selected = "None"
      )
      
    )
    
  })
  
  
  #############################################################
  # UPDATE COVARIATE SELECTION
  #############################################################
  
  observe({
    
    df <- data_input()
    
    updateSelectizeInput(
      session,
      "covariates",
      choices = setdiff(names(df), c(input$time_var, input$event_var)),
      server = TRUE
    )
    
  })
  
  
  #############################################################
  # SURVIVAL OBJECT
  #############################################################
  
  # Build survival model using standardized column names
  surv_object <- reactive({
    
    df <- data_input()
    
    req(input$time_var, input$event_var)
    
    # create standardized variables to avoid NSE problems
    df$.time  <- df[[input$time_var]]
    df$.event <- df[[input$event_var]]
    
    if (input$group_var == "None") {
      
      fit <- survfit(Surv(.time, .event) ~ 1, data = df)
      
    } else {
      
      df$.group <- as.factor(df[[input$group_var]])
      fit <- survfit(Surv(.time, .event) ~ .group, data = df)
      
    }
    
    list(fit = fit, data = df)
    
  })
  
  
  #############################################################
  # KAPLAN-MEIER CURVE
  #############################################################
  
  output$km_plot <- renderPlot({
    
    obj <- surv_object()
    
    p <- ggsurvplot(
      obj$fit,
      data = obj$data,
      conf.int = TRUE,
      pval = input$group_var != "None",
      risk.table = FALSE,
      xlab = "Time",
      ylab = "Survival Probability",
      ggtheme = theme_minimal()
    )
    
    print(p)
    
  })
  
  
  #############################################################
  # DESCRIPTIVE SURVIVAL STATISTICS
  #############################################################
  
  output$surv_stats <- renderTable({
    
    obj <- surv_object()
    fit <- obj$fit
    df  <- obj$data
    
    fit_table <- summary(fit)$table
    
    # Case: no stratification
    if (is.null(dim(fit_table))) {
      
      total_n <- nrow(df)
      total_events <- sum(df$.event)
      
      med <- as.numeric(fit_table["median"])
      surv_365 <- summary(fit, times = 365, extend = TRUE)$surv[1]
      
      data.frame(
        Group = "Overall",
        Total_Patients = total_n,
        Events = total_events,
        Event_Rate_Percent = round(100 * total_events / total_n,1),
        Median_Survival_Time = round(med,1),
        Survival_Probability_365_Days = round(surv_365,3)
      )
      
    } else {
      
      med <- as.numeric(fit_table[, "median"])
      surv_365 <- summary(fit, times = 365, extend = TRUE)$surv
      
      data.frame(
        Group = rownames(fit_table),
        Median_Survival_Time = round(med,1),
        Survival_Probability_365_Days = round(surv_365,3)
      )
      
    }
    
  })
  
  
  #############################################################
  # UNIVARIATE COX MODEL
  #############################################################
  
  output$cox_results <- renderTable({
    
    if (input$group_var == "None") {
      return(data.frame(
        Message = "Select a grouping variable to run Cox regression."
      ))
    }
    
    obj <- surv_object()
    df  <- obj$data
    
    model <- coxph(Surv(.time, .event) ~ .group, data = df)
    
    s <- summary(model)
    
    data.frame(
      Variable = rownames(s$coefficients),
      Hazard_Ratio = round(s$coefficients[,"exp(coef)"],3),
      CI_95 = paste0(
        round(s$conf.int[,"lower .95"],3),
        " - ",
        round(s$conf.int[,"upper .95"],3)
      ),
      P_Value = signif(s$coefficients[,"Pr(>|z|)"],3)
    )
    
  })
  
  
  #############################################################
  # COX MODEL INTERPRETATION
  #############################################################
  
  output$cox_interpretation <- renderPrint({
    
    if (input$group_var == "None") return(NULL)
    
    obj <- surv_object()
    df  <- obj$data
    
    model <- coxph(Surv(.time, .event) ~ .group, data = df)
    
    s <- summary(model)
    
    hr <- s$coefficients[,"exp(coef)"]
    p  <- s$coefficients[,"Pr(>|z|)"]
    
    variable <- gsub(".group","",rownames(s$coefficients), fixed = TRUE)
    
    significance <- ifelse(
      p < 0.05,
      "statistically significant",
      "not statistically significant"
    )
    
    cat(
      "Interpretation:\n\n",
      variable,
      "has a hazard ratio of",
      round(hr,3),
      "relative to the reference group.\n",
      "This effect is",
      significance,
      "(p =", signif(p,3), ")."
    )
    
  })
  
  
  #############################################################
  # MULTIVARIABLE COX MODEL
  #############################################################
  
  output$cox_multi_results <- renderTable({
    
    if (input$group_var == "None") return(NULL)
    if (length(input$covariates) == 0) return(NULL)
    
    obj <- surv_object()
    df  <- obj$data
    
    formula_str <- paste(
      "Surv(.time, .event) ~ .group +",
      paste(input$covariates, collapse = " + ")
    )
    
    model <- coxph(as.formula(formula_str), data = df)
    
    s <- summary(model)
    
    data.frame(
      Variable = rownames(s$coefficients),
      Hazard_Ratio = round(s$coefficients[,"exp(coef)"],3),
      CI_95 = paste0(
        round(s$conf.int[,"lower .95"],3),
        " - ",
        round(s$conf.int[,"upper .95"],3)
      ),
      P_Value = signif(s$coefficients[,"Pr(>|z|)"],3)
    )
    
  })
  
  
  #############################################################
  # FOREST PLOT
  #############################################################
  
  output$forest_plot <- renderPlot({
    
    if (input$group_var == "None") return(NULL)
    if (length(input$covariates) == 0) return(NULL)
    
    obj <- surv_object()
    df  <- obj$data
    
    formula_str <- paste(
      "Surv(.time, .event) ~ .group +",
      paste(input$covariates, collapse = " + ")
    )
    
    model <- coxph(as.formula(formula_str), data = df)
    
    ggforest(model, data = df)
    
  })
  
  
  #############################################################
  # PROPORTIONAL HAZARDS TEST
  #############################################################
  
  output$ph_test <- renderTable({
    
    if (input$group_var == "None") return(NULL)
    if (length(input$covariates) == 0) return(NULL)
    
    obj <- surv_object()
    df  <- obj$data
    
    formula_str <- paste(
      "Surv(.time, .event) ~ .group +",
      paste(input$covariates, collapse = " + ")
    )
    
    model <- coxph(as.formula(formula_str), data = df)
    
    ph <- cox.zph(model)
    
    data.frame(
      Variable = rownames(ph$table),
      Chi_Square = round(ph$table[,"chisq"],3),
      P_Value = signif(ph$table[,"p"],3)
    )
    
  })
  
  
  #############################################################
  # SCHOENFELD RESIDUAL PLOTS
  #############################################################
  
  output$ph_plot <- renderPlot({
    
    if (input$group_var == "None") return(NULL)
    if (length(input$covariates) == 0) return(NULL)
    
    obj <- surv_object()
    df  <- obj$data
    
    formula_str <- paste(
      "Surv(.time, .event) ~ .group +",
      paste(input$covariates, collapse = " + ")
    )
    
    model <- coxph(as.formula(formula_str), data = df)
    
    ph <- cox.zph(model)
    
    plot(ph)
    
  })
  
  
  #############################################################
  # PH INTERPRETATION
  #############################################################
  
  output$ph_interpretation <- renderPrint({
    
    if (input$group_var == "None") return(NULL)
    if (length(input$covariates) == 0) return(NULL)
    
    obj <- surv_object()
    df  <- obj$data
    
    formula_str <- paste(
      "Surv(.time, .event) ~ .group +",
      paste(input$covariates, collapse = " + ")
    )
    
    model <- coxph(as.formula(formula_str), data = df)
    
    ph <- cox.zph(model)
    
    global_p <- ph$table[nrow(ph$table),"p"]
    
    if(global_p < 0.05){
      
      cat(
        "Global PH test p-value:",
        signif(global_p,3),
        "\nThe proportional hazards assumption may be violated."
      )
      
    } else {
      
      cat(
        "Global PH test p-value:",
        signif(global_p,3),
        "\nNo evidence against the proportional hazards assumption."
      )
      
    }
    
  })
  
}

###############################################################
# RUN APP
###############################################################

shinyApp(ui = ui, server = server)