# Interactive Survival Explorer

A Shiny application for interactive survival analysis using Kaplan–Meier estimation and Cox proportional hazards modeling.

This project demonstrates how statistical survival models can be explored interactively using **R**, **Shiny**, and the **survival** ecosystem. It is designed as a portfolio project for **bioinformatics and biostatistics workflows**.

---

## Overview

Survival analysis is widely used in biomedical research to model time-to-event outcomes such as patient survival, disease recurrence, or treatment response.

This application allows users to:

* Explore **Kaplan–Meier survival curves**
* Perform **Cox proportional hazards regression**
* Fit **multivariable survival models**
* Visualize **hazard ratios with forest plots**
* Test the **proportional hazards assumption**
* Upload custom datasets for interactive analysis

The app ships with the `lung` dataset from the **survival** package as an example dataset.

---

## Features

### Kaplan–Meier Survival Analysis

* Survival curves with confidence intervals
* Optional stratification by group
* Log-rank test p-value

### Descriptive Survival Statistics

* Total number of patients
* Number of events
* Event rate
* Median survival time
* Survival probability at 365 days

### Cox Proportional Hazards Models

* Univariate Cox regression
* Hazard ratios with 95% confidence intervals
* Automatic interpretation text

### Multivariable Cox Models

* User-selectable covariates
* Forest plot visualization of hazard ratios
* Adjusted hazard ratios for multiple predictors

### Proportional Hazards Diagnostics

* `cox.zph()` proportional hazards test
* Schoenfeld residual plots
* Global PH assumption interpretation

### Dynamic Dataset Support

* Upload custom CSV datasets
* Automatically select survival variables
* Run survival models interactively

---

## Example Dataset

The application uses the `lung` dataset from the **survival** package by default.

This dataset contains clinical data from lung cancer patients, including:

* survival time
* censoring status
* sex
* age
* performance scores
* weight loss
* calorie intake

---

## Technologies Used

* **R**
* **Shiny**
* **survival**
* **survminer**
* **dplyr**

These packages provide tools for statistical modeling, visualization, and interactive dashboards.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/interactive-survival-explorer.git
cd interactive-survival-explorer
```

Install required packages in R:

```r
install.packages(c(
  "shiny",
  "survival",
  "survminer",
  "dplyr"
))
```

Run the app:

```r
shiny::runApp()
```

---

## Project Structure

```
interactive-survival-explorer/
│
├── app.R
└── README.md

```

* **app.R** – main Shiny application
* **README.md** – project documentation

---

## Example Workflow

1. Launch the Shiny app
2. Select time and event variables
3. Choose a grouping variable for stratified survival analysis
4. Select covariates for multivariable Cox regression
5. Examine survival curves, hazard ratios, and PH diagnostics

---

## Applications in Bioinformatics

Survival analysis is frequently used in:

* cancer genomics
* biomarker discovery
* clinical trial analysis
* treatment response modeling
* patient risk stratification

This project demonstrates how these statistical tools can be deployed in **interactive analytical environments**.

---

## Future Improvements

Potential extensions for this project include:

* gene expression survival analysis
* TCGA-style high/low expression survival comparisons
* LASSO Cox regression for feature selection
* downloadable survival reports
* improved dashboard UI
* optional example datasets

---

## Author

Sergio Fernández Álvarez
Bioinformatics / Biostatistics portfolio project.
