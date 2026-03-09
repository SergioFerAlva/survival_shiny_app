###############################################################
# Interactive Survival Explorer
# Biostatistics / Bioinformatics Portfolio App
###############################################################

library(shiny)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(bslib)

###############################################################
# EXAMPLE DATASETS
###############################################################

lung <- survival::lung %>%
  mutate(
    status = as.numeric(status == 2),
    sex = factor(sex, labels = c("Male","Female"))
  )

pbc <- survival::pbc %>%
  mutate(
    status = as.numeric(status == 2)
  )

veteran <- survival::veteran %>%
  mutate(
    status = as.numeric(status == 1)
  )

###############################################################
# USER INTERFACE
###############################################################

ui <- fluidPage(
  
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  titlePanel("Interactive Survival Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h4("Dataset"),
      
      selectInput(
        "dataset",
        "Choose example dataset",
        choices = c("Lung","PBC","Veteran")
      ),
      
      fileInput(
        "file",
        "Upload CSV dataset",
        accept = ".csv"
      ),
      
      hr(),
      
      h4("Variables"),
      
      uiOutput("var_select"),
      
      selectizeInput(
        "covariates",
        "Covariates for multivariable Cox model",
        choices = NULL,
        multiple = TRUE
      ),
      
      hr(),
      
      h4("Plot Options"),
      
      checkboxInput("show_ci","Show confidence interval", TRUE),
      checkboxInput("risk_table","Show risk table", TRUE),
      
      hr(),
      
      downloadButton("download_cox","Download Cox Results"),
      downloadButton("download_km","Download KM Plot")
      
    ),
    
    mainPanel(
      
      tabsetPanel(
        
        tabPanel(
          "Summary",
          
          h3("Dataset Summary"),
          tableOutput("dataset_summary"),
          
          h3("Descriptive Survival Statistics"),
          tableOutput("surv_stats"),
          
          h3("Log-Rank Test"),
          tableOutput("logrank_test")
          
        ),
        
        tabPanel(
          "Kaplan-Meier",
          
          plotOutput("km_plot")
          
        ),
        
        tabPanel(
          "Cox Model",
          
          h3("Univariate Cox Model"),
          tableOutput("cox_results"),
          verbatimTextOutput("cox_interpretation"),
          
          h3("Multivariable Cox Model"),
          tableOutput("cox_multi_results"),
          
          h3("Forest Plot"),
          plotOutput("forest_plot"),
          
          h3("Model Performance"),
          tableOutput("model_performance")
          
        ),
        
        tabPanel(
          "Diagnostics",
          
          h3("Proportional Hazards Test"),
          tableOutput("ph_test"),
          
          plotOutput("ph_plot"),
          
          verbatimTextOutput("ph_interpretation")
          
        ),
        
        tabPanel(
          "Methods",
          
          h3("Statistical Methods"),
          
          tags$ul(
            tags$li("Kaplan-Meier survival estimation"),
            tags$li("Log-rank test for group comparison"),
            tags$li("Cox proportional hazards regression"),
            tags$li("Proportional hazards assumption testing")
          )
          
        )
      )
    )
  )
)

###############################################################
# SERVER
###############################################################

server <- function(input, output, session) {
  
  ###############################################################
  # DATA INPUT
  ###############################################################
  
  data_input <- reactive({
    
    if(!is.null(input$file)){
      return(read.csv(input$file$datapath))
    }
    
    if(input$dataset == "Lung") return(lung)
    if(input$dataset == "PBC") return(pbc)
    if(input$dataset == "Veteran") return(veteran)
    
  })
  
  ###############################################################
  # VARIABLE SELECTION
  ###############################################################
  
  output$var_select <- renderUI({
    
    df <- data_input()
    
    time_default <- if("time" %in% names(df)) "time" else names(df)[1]
    event_default <- if("status" %in% names(df)) "status" else names(df)[2]
    
    tagList(
      
      selectInput(
        "time_var",
        "Time variable:",
        choices = names(df),
        selected = time_default
      ),
      
      selectInput(
        "event_var",
        "Event variable:",
        choices = names(df),
        selected = event_default
      ),
      
      selectInput(
        "group_var",
        "Grouping variable:",
        choices = c("None", names(df)),
        selected = "None"
      )
      
    )
    
  })
  
  ###############################################################
  # UPDATE COVARIATES
  ###############################################################
  
  observe({
    
    df <- data_input()
    
    covars <- setdiff(names(df), c(input$time_var, input$event_var))
    
    updateSelectizeInput(
      session,
      "covariates",
      choices = covars,
      server = TRUE
    )
    
  })
  
  ###############################################################
  # SURVIVAL OBJECT
  ###############################################################
  
  surv_object <- reactive({
    
    df <- data_input()
    
    df$.time <- df[[input$time_var]]
    df$.event <- df[[input$event_var]]
    
    if(input$group_var != "None"){
      
      df$.group <- as.factor(df[[input$group_var]])
      fit <- survfit(Surv(.time,.event)~.group,data=df)
      
    } else {
      
      fit <- survfit(Surv(.time,.event)~1,data=df)
      
    }
    
    list(fit=fit,data=df)
    
  })
  
  ###############################################################
  # DATASET SUMMARY
  ###############################################################
  
  output$dataset_summary <- renderTable({
    
    df <- data_input()
    
    data.frame(
      Total_Patients = nrow(df),
      Variables = ncol(df),
      Missing_Values = sum(is.na(df))
    )
    
  })
  
  ###############################################################
  # SURVIVAL STATISTICS (GROUP SAFE)
  ###############################################################
  
  output$surv_stats <- renderTable({
    
    obj <- surv_object()
    fit <- obj$fit
    
    med <- summary(fit)$table
    
    surv_times <- summary(fit, times=c(365,1095,1825), extend=TRUE)
    
    if(is.null(dim(med))){
      
      data.frame(
        Group="Overall",
        Median_Survival=round(med["median"],1),
        Survival_1yr=round(surv_times$surv[1],3),
        Survival_3yr=round(surv_times$surv[2],3),
        Survival_5yr=round(surv_times$surv[3],3)
      )
      
    } else {
      
      data.frame(
        Group=rownames(med),
        Median_Survival=round(med[,"median"],1)
      )
      
    }
    
  })
  
  ###############################################################
  # LOG-RANK TEST
  ###############################################################
  
  output$logrank_test <- renderTable({
    
    if(input$group_var=="None") return(NULL)
    
    df <- data_input()
    
    df$.time <- df[[input$time_var]]
    df$.event <- df[[input$event_var]]
    df$.group <- as.factor(df[[input$group_var]])
    
    lr <- survdiff(Surv(.time,.event)~.group,data=df)
    
    data.frame(
      Chi_Square = round(lr$chisq,3),
      P_Value = signif(pchisq(lr$chisq,length(lr$n)-1,lower.tail=FALSE),3)
    )
    
  })
  
  ###############################################################
  # KM PLOT
  ###############################################################
  
  output$km_plot <- renderPlot({
    
    obj <- surv_object()
    
    p <- ggsurvplot(
      obj$fit,
      data=obj$data,
      conf.int=input$show_ci,
      risk.table=input$risk_table,
      ggtheme=theme_minimal()
    )
    
    print(p)
    
  })
  
  ###############################################################
  # COX MODEL
  ###############################################################
  
  output$cox_results <- renderTable({
    
    if(input$group_var=="None") return(NULL)
    
    obj <- surv_object()
    
    model <- coxph(Surv(.time,.event)~.group,data=obj$data)
    
    s <- summary(model)
    
    data.frame(
      Variable=rownames(s$coefficients),
      Hazard_Ratio=round(s$coefficients[,"exp(coef)"],3),
      CI_Lower=round(s$conf.int[,"lower .95"],3),
      CI_Upper=round(s$conf.int[,"upper .95"],3),
      P_Value=signif(s$coefficients[,"Pr(>|z|)"],3)
    )
    
  })
  
  ###############################################################
  # COX INTERPRETATION
  ###############################################################
  
  output$cox_interpretation <- renderPrint({
    
    if(input$group_var=="None") return(NULL)
    
    obj <- surv_object()
    
    model <- coxph(Surv(.time,.event)~.group,data=obj$data)
    
    s <- summary(model)
    
    hr <- s$coefficients[,"exp(coef)"]
    p <- s$coefficients[,"Pr(>|z|)"]
    
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
      "relative to the reference group.\n\n",
      "This effect is",
      significance,
      "(p =", signif(p,3), ")."
    )
    
  })
  
  ###############################################################
  # MULTIVARIABLE COX
  ###############################################################
  
  output$cox_multi_results <- renderTable({
    
    if(length(input$covariates)==0) return(NULL)
    
    obj <- surv_object()
    
    formula_str <- paste(
      "Surv(.time,.event)~",
      paste(input$covariates,collapse="+")
    )
    
    model <- coxph(as.formula(formula_str),data=obj$data)
    
    s <- summary(model)
    
    data.frame(
      Variable=rownames(s$coefficients),
      Hazard_Ratio=round(s$coefficients[,"exp(coef)"],3),
      P_Value=signif(s$coefficients[,"Pr(>|z|)"],3)
    )
    
  })
  
  ###############################################################
  # FOREST PLOT
  ###############################################################
  
  output$forest_plot <- renderPlot({
    
    if(length(input$covariates)==0) return(NULL)
    
    obj <- surv_object()
    
    formula_str <- paste(
      "Surv(.time,.event)~",
      paste(input$covariates,collapse="+")
    )
    
    model <- coxph(as.formula(formula_str),data=obj$data)
    
    ggforest(model,data=obj$data)
    
  })
  
  ###############################################################
  # MODEL PERFORMANCE
  ###############################################################
  
  output$model_performance <- renderTable({
    
    if(length(input$covariates)==0) return(NULL)
    
    obj <- surv_object()
    
    formula_str <- paste(
      "Surv(.time,.event)~",
      paste(input$covariates,collapse="+")
    )
    
    model <- coxph(as.formula(formula_str),data=obj$data)
    
    s <- summary(model)
    
    data.frame(
      C_Index = round(s$concordance[1],3),
      Likelihood_Ratio_p = signif(s$logtest["pvalue"],3)
    )
    
  })
  
  ###############################################################
  # PH TEST
  ###############################################################
  
  output$ph_test <- renderTable({
    
    if(length(input$covariates)==0) return(NULL)
    
    obj <- surv_object()
    
    formula_str <- paste(
      "Surv(.time,.event)~",
      paste(input$covariates,collapse="+")
    )
    
    model <- coxph(as.formula(formula_str),data=obj$data)
    
    ph <- cox.zph(model)
    
    data.frame(
      Variable=rownames(ph$table),
      P_Value=signif(ph$table[,"p"],3)
    )
    
  })
  
  ###############################################################
  # PH PLOT
  ###############################################################
  
  output$ph_plot <- renderPlot({
    
    if(length(input$covariates)==0) return(NULL)
    
    obj <- surv_object()
    
    formula_str <- paste(
      "Surv(.time,.event)~",
      paste(input$covariates,collapse="+")
    )
    
    model <- coxph(as.formula(formula_str),data=obj$data)
    
    ph <- cox.zph(model)
    
    plot(ph)
    
  })
  
  ###############################################################
  # PH INTERPRETATION
  ###############################################################
  
  output$ph_interpretation <- renderPrint({
    
    if(length(input$covariates)==0) return(NULL)
    
    obj <- surv_object()
    
    formula_str <- paste(
      "Surv(.time,.event)~",
      paste(input$covariates,collapse="+")
    )
    
    model <- coxph(as.formula(formula_str),data=obj$data)
    
    ph <- cox.zph(model)
    
    p <- ph$table[nrow(ph$table),"p"]
    
    if(p < 0.05){
      cat("Possible violation of proportional hazards assumption.")
    } else {
      cat("No evidence of PH violation.")
    }
    
  })
  
}

###############################################################
# RUN APP
###############################################################

shinyApp(ui=ui,server=server)