# Interactive Survival Explorer

An interactive **Shiny application for survival analysis** using Kaplan–Meier estimation and Cox proportional hazards modeling.

This project demonstrates how statistical survival methods commonly used in **biostatistics, clinical research, and bioinformatics** can be explored through an interactive dashboard built with **R**, **Shiny**, and the **survival analysis ecosystem**.

The goal of this project is to showcase applied statistical modeling, reproducible analysis workflows, and interactive scientific tools suitable for **biostatistics and bioinformatics portfolios**.

---

# Overview

Survival analysis is widely used in biomedical research to model **time-to-event outcomes**, such as:

* patient survival
* disease recurrence
* treatment response
* time to relapse

This application allows users to perform common survival analyses interactively without writing code.

Users can:

* explore Kaplan–Meier survival curves
* compare groups using the log-rank test
* fit Cox proportional hazards models
* visualize hazard ratios with forest plots
* evaluate model assumptions
* upload custom datasets for analysis

---

# Features

## Kaplan–Meier Survival Analysis

* Kaplan–Meier survival curves
* Optional confidence intervals
* Risk table visualization
* Group stratification
* Log-rank test for survival differences

---

## Descriptive Survival Statistics

The app provides summary statistics including:

* number of patients
* number of events
* median survival time
* survival probabilities at:

  * 1 year
  * 3 years
  * 5 years

---

## Cox Proportional Hazards Models

The application supports both **univariate and multivariable Cox regression**.

Features include:

* hazard ratios (HR)
* 95% confidence intervals
* p-values
* automatic interpretation of hazard ratios

---

## Multivariable Cox Models

Users can select multiple covariates to fit adjusted survival models.

The app provides:

* adjusted hazard ratios
* forest plot visualization
* model performance metrics including:

  * concordance index (C-index)
  * likelihood ratio test

---

## Model Diagnostics

The application includes diagnostics for evaluating model assumptions.

These include:

* proportional hazards test (`cox.zph`)
* Schoenfeld residual plots
* global interpretation of PH assumption validity

---

## Interactive Visualization Options

Users can customize survival plots by:

* enabling/disabling confidence intervals
* toggling risk tables
* stratifying survival curves by groups

---

## Dynamic Dataset Support

The application supports both **built-in example datasets** and **user-uploaded data**.

Example datasets included:

* `lung` — lung cancer survival dataset
* `pbc` — primary biliary cirrhosis dataset
* `veteran` — veteran lung cancer trial dataset

Users can also upload **custom CSV datasets** and select:

* time variable
* event indicator
* grouping variables
* model covariates

---

# Technologies Used

This project is built using the following R ecosystem tools:

* **R**
* **Shiny**
* **survival**
* **survminer**
* **dplyr**
* **ggplot2**
* **bslib**

These packages provide statistical modeling, visualization, and interactive web application capabilities.

---

# Installation

Clone the repository:

```
git clone https://github.com/yourusername/interactive-survival-explorer.git
cd interactive-survival-explorer
```

Install required packages in R:

```r
install.packages(c(
  "shiny",
  "survival",
  "survminer",
  "dplyr",
  "ggplot2",
  "bslib"
))
```

Run the application:

```r
shiny::runApp()
```

---

# Project Structure

```
interactive-survival-explorer/
│
├── app.R
├── README.md
├── .gitignore
```

* **app.R** — main Shiny application
* **README.md** — project documentation

---

# Example Workflow

A typical workflow in the application:

1. Launch the Shiny app
2. Select an example dataset or upload a CSV dataset
3. Choose the survival **time** and **event** variables
4. Select a **grouping variable** for Kaplan–Meier comparison
5. Run **Cox proportional hazards regression**
6. Evaluate model performance and diagnostics

---

# Applications in Bioinformatics and Clinical Research

Survival analysis plays a critical role in many areas of biomedical data science, including:

* cancer genomics
* biomarker discovery
* clinical trial analysis
* treatment response modeling
* patient risk stratification
* translational bioinformatics

This project demonstrates how these statistical tools can be deployed in **interactive analytical environments for exploratory research and reproducible workflows**.

---

# Future Improvements

Potential future extensions include:

* TCGA gene expression survival analysis
* high vs low gene expression Kaplan–Meier analysis
* penalized Cox models (LASSO survival)
* time-dependent covariates
* downloadable survival analysis reports

---

# Author

**Sergio Fernández Álvarez**

Bioinformatics / Biostatistics portfolio project.
