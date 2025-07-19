#
# Shiny App for Population Bioequivalence (PBE) Test
# Based on the FDA Guidance: "Statistical Approaches to Establishing Bioequivalence" (December 2022)
# Version 5: Added a detailed, step-by-step calculation breakdown with formulas.
#

# Load the necessary libraries
# 'shiny' is for the web app framework
# 'readxl' is for reading Excel files (.xls, .xlsx)
# 'dplyr' is used for easier data manipulation
library(shiny)
library(readxl)
library(dplyr)

# Define the user interface (UI)
ui <- fluidPage(
  # Add MathJax support to render LaTeX formulas
  withMathJax(),
  titlePanel("Population Bioequivalence (PBE) Test Tool"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Instructions"),
      
      # 1. File Upload
      fileInput("file1", "1. Upload Your Data File",
                multiple = FALSE,
                accept = c(".csv", ".xls", ".xlsx",
                           "text/csv",
                           "text/comma-separated-values,text/plain",
                           "application/vnd.ms-excel",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")),
      
      # 2. Data Format Instructions
      h5("Required Data Format:"),
      p("The app supports long-format data. After uploading, options will appear below to map your columns (e.g., select your 'Outcome' column). Your 'Outcome' column should contain log-transformed data."),
      
      # 3. Dynamic UI for column selection will appear here
      uiOutput("data_mapping_ui"),
      
      # 4. Action Button
      tags$hr(),
      # The button is disabled until a file is loaded and columns are mapped
      uiOutput("run_button_ui"),
      
      # 5. Citation
      tags$hr(),
      p(strong("Source:"), "U.S. FDA, Guidance for Industry, 'Statistical Approaches to Establishing Bioequivalence', Rev. 1, Dec 2022.")
    ),
    
    mainPanel(
      # UI to display results
      uiOutput("results_ui")
    )
  )
)
# Define the server logic
server <- function(input, output, session) {
  
  # Reactive expression to read the uploaded file
  raw_data <- reactive({
    req(input$file1)
    
    tryCatch({
      ext <- tools::file_ext(input$file1$name)
      if (ext == "csv") {
        read.csv(input$file1$datapath, header = TRUE, stringsAsFactors = FALSE)
      } else if (ext %in% c("xls", "xlsx")) {
        readxl::read_excel(input$file1$datapath, col_names = TRUE)
      } else {
        stop("Unsupported file type. Please upload a CSV or Excel file.")
      }
    }, error = function(e) {
      shiny::validate(need(FALSE, paste("File Reading Error:", e$message)))
    })
  })
  
  # Render dynamic UI for data mapping based on the uploaded file
  output$data_mapping_ui <- renderUI({
    df <- raw_data()
    if (is.null(df)) return(NULL)
    
    col_names <- names(df)
    
    tagList(
      tags$hr(),
      h5("3. Map Your Data Columns"),
      p("Select the columns from your file that correspond to the required data."),
      
      selectInput("treatment_col", "Select Treatment Group Column:", choices = col_names, selected = col_names[1]),
      
      p("The 'Outcome' is the measured result (e.g., log(AUC) or log(Cmax))."),
      selectInput("outcome_col", "Select Outcome (Value) Column:", choices = col_names, selected = col_names[2]),
      
      uiOutput("treatment_level_selectors_ui")
    )
  })
  
  # Render selectors for Test/Reference levels based on the chosen treatment column
  output$treatment_level_selectors_ui <- renderUI({
    df <- raw_data()
    req(input$treatment_col)
    
    treatment_levels <- unique(df[[input$treatment_col]])
    
    tagList(
      selectInput("test_level", "Which value is 'Test'?", choices = treatment_levels, selected = treatment_levels[1]),
      selectInput("ref_level", "Which value is 'Reference'?", choices = treatment_levels, selected = treatment_levels[2])
    )
  })
  
  # Render the run button dynamically
  output$run_button_ui <- renderUI({
    if (is.null(raw_data())) return(NULL)
    actionButton("run_test", "Run PBE Test", class = "btn-primary")
  })
  
  # Reactive value to store the results
  results <- eventReactive(input$run_test, {
    
    req(raw_data(), input$treatment_col, input$outcome_col, input$test_level, input$ref_level)
    validate(
      need(input$test_level != input$ref_level, "Error: 'Test' and 'Reference' product levels cannot be the same.")
    )
    
    df <- raw_data()
    T_data <- df[[input$outcome_col]][df[[input$treatment_col]] == input$test_level]
    R_data <- df[[input$outcome_col]][df[[input$treatment_col]] == input$ref_level]
    
    if (length(T_data) == 0 || length(R_data) == 0) {
      stop("Error: No data found for the selected 'Test' or 'Reference' levels.")
    }
    if (!is.numeric(T_data) || !is.numeric(R_data)) {
      stop("Error: The selected 'Outcome' column must contain numeric data.")
    }
    
    T_data_clean <- T_data[!is.na(T_data)]
    R_data_clean <- R_data[!is.na(R_data)]
    
    theta_P <- 2.089
    sigma0_sq <- 0.01
    alpha <- 0.05
    
    n_T <- length(T_data_clean)
    n_R <- length(R_data_clean)
    
    mean_T <- mean(T_data_clean)
    mean_R <- mean(R_data_clean)
    
    var_T <- var(T_data_clean)
    var_R <- var(R_data_clean)
    
    gamma1_hat <- (mean_T - mean_R)^2
    gamma2_hat <- var_T
    
    if (var_R > sigma0_sq) {
      gamma3_hat <- var_R + theta_P * var_R
    } else {
      gamma3_hat <- var_R + theta_P * sigma0_sq
    }
    
    se_diff <- sqrt(var_T / n_T + var_R / n_R)
    df_welch <- (var_T / n_T + var_R / n_R)^2 / (((var_T / n_T)^2 / (n_T - 1)) + ((var_R / n_R)^2 / (n_R - 1)))
    t_crit <- qt(1 - alpha, df_welch)
    ci_diff <- (mean_T - mean_R) + c(-1, 1) * t_crit * se_diff
    gamma1_tilde <- (max(abs(ci_diff)))^2
    
    gamma2_tilde <- (n_T - 1) * var_T / qchisq(alpha, df = n_T - 1)
    
    lower_bound_var_R <- (n_R - 1) * var_R / qchisq(1 - alpha, df = n_R - 1)
    if (var_R > sigma0_sq) {
      gamma3_tilde <- lower_bound_var_R + theta_P * lower_bound_var_R
    } else {
      gamma3_tilde <- lower_bound_var_R + theta_P * sigma0_sq
    }
    
    gamma_U_hat <- (gamma1_hat + gamma2_hat - gamma3_hat) + 
      sqrt((gamma1_tilde - gamma1_hat)^2 + 
             (gamma2_tilde - gamma2_hat)^2 + 
             (gamma3_tilde - gamma3_hat)^2)
    
    pbe_passed <- gamma_U_hat <= 0
    
    list(
      summary_stats = data.frame(
        Parameter = c("N (Count)", "Mean", "Variance"),
        Test = c(n_T, round(mean_T, 5), round(var_T, 5)),
        Reference = c(n_R, round(mean_R, 5), round(var_R, 5))
      ),
      gamma_U_hat = gamma_U_hat,
      pbe_passed = pbe_passed,
      calculation_steps = list(
        var_R_val = var_R,
        gamma1_hat = gamma1_hat,
        gamma2_hat = gamma2_hat,
        gamma3_hat = gamma3_hat,
        gamma1_tilde = gamma1_tilde,
        gamma2_tilde = gamma2_tilde,
        gamma3_tilde = gamma3_tilde
      )
    )
  })
  
  # Render the results in the main panel
  output$results_ui <- renderUI({
    req(input$run_test)
    
    res <- results()
    if (is.null(res)) return()
    
    result_text <- if (res$pbe_passed) "SUPPORTED" else "NOT SUPPORTED"
    result_color <- if (res$pbe_passed) "color:green;" else "color:red;"
    
    tagList(
      h3("PBE Test Results"),
      tags$hr(),
      h4("1. Summary of Processed Data"),
      p(paste0("(Using ", res$summary_stats$Test[1], " Test and ", res$summary_stats$Reference[1], " Reference observations after removing missing values)")),
      tableOutput("summary_table"),
      tags$hr(),
      h4("2. Detailed Calculation Steps"),
      p("The following steps detail the calculation of the linearized criterion (γ) and its 95% upper confidence bound, based on the FDA guidance."),
      uiOutput("calculation_details"),
      tags$hr(),
      h4("3. PBE Test Criterion"),
      p("The test requires the 95% upper confidence bound of the linearized criterion (γ) to be less than or equal to zero."),
      div(
        style = "padding: 15px; border: 1px solid #ccc; border-radius: 5px; background-color: #f9f9f9; text-align: center;",
        h5("Calculated 95% Upper Bound (γ̂_U):"),
        h3(strong(round(res$gamma_U_hat, 5))),
        h5("Criterion:"),
        h3(strong("γ̂_U ≤ 0"))
      ),
      tags$hr(),
      h4("4. Final Conclusion"),
      div(
        style = paste("padding: 20px; border: 2px solid; text-align: center; font-size: 24px; font-weight: bold;", result_color),
        p("Population Bioequivalence is ", strong(result_text))
      )
    )
  })
  
  # Render the summary table
  output$summary_table <- renderTable({
    req(input$run_test)
    results()$summary_stats
  })
  
  # *** CORRECTED SECTION ***
  # Render the detailed calculation steps with correctly formatted MathJax
  output$calculation_details <- renderUI({
    req(input$run_test)
    res <- results()
    steps <- res$calculation_steps
    var_R_val <- steps$var_R_val
    sigma0_sq <- 0.01
    
    # Define the conditional text and formulas
    condition_text <- paste0("(Since σ̂²_R of ", round(var_R_val, 4), 
                             (if(var_R_val > sigma0_sq) " > " else " ≤ "), 
                             sigma0_sq, ", the scaled formula is used.)")
    
    gamma3_formula <- if (var_R_val > sigma0_sq) {
      "$$\\hat{\\gamma}_3 = \\hat{\\sigma}_R^2 + \\theta_P \\cdot \\hat{\\sigma}_R^2$$"
    } else {
      "$$\\hat{\\gamma}_3 = \\hat{\\sigma}_R^2 + \\theta_P \\cdot \\sigma_0^2$$"
    }
    
    gamma3_tilde_formula <- if (var_R_val > sigma0_sq) {
      "$$\\tilde{\\gamma}_3 = \\text{LowerCI}(\\sigma_R^2) + \\theta_P \\cdot \\text{LowerCI}(\\sigma_R^2)$$"
    } else {
      "$$\\tilde{\\gamma}_3 = \\text{LowerCI}(\\sigma_R^2) + \\theta_P \\cdot \\sigma_0^2$$"
    }
    
    # Wrap the entire output in withMathJax() to ensure it gets rendered
    withMathJax(
      div(
        style="padding: 10px; border-left: 3px solid #eee;",
        h5(strong("Part A: Point Estimators (hat values)")),
        p(strong("1. γ̂₁: "), sprintf("$$\\hat{\\gamma}_1 = (\\hat{\\mu}_T - \\hat{\\mu}_R)^2 = %.5f$$", steps$gamma1_hat)),
        p(strong("2. γ̂₂: "), sprintf("$$\\hat{\\gamma}_2 = \\hat{\\sigma}_T^2 = %.5f$$", steps$gamma2_hat)),
        p(strong("3. γ̂₃: "), condition_text, sprintf("%s = %.5f", gamma3_formula, steps$gamma3_hat)),
        
        tags$hr(),
        h5(strong("Part B: Confidence Bounds (tilde values)")),
        p(strong("4. γ̃₁ (95% Upper): "), sprintf("$$\\tilde{\\gamma}_1 = (\\text{max}|\\text{CI}(\\mu_T - \\mu_R)|)^2 = %.5f$$", steps$gamma1_tilde)),
        p(strong("5. γ̃₂ (95% Upper): "), sprintf("$$\\tilde{\\gamma}_2 = \\frac{(n_T-1)\\hat{\\sigma}_T^2}{\\chi^2_{\\alpha, n_T-1}} = %.5f$$", steps$gamma2_tilde)),
        p(strong("6. γ̃₃ (95% Lower): "), sprintf("%s = %.5f", gamma3_tilde_formula, steps$gamma3_tilde)),
        
        tags$hr(),
        h5(strong("Part C: Final Calculation of the Upper Bound (γ̂_U)")),
        p("This combines the values from Parts A and B:"),
        p("$$\\hat{\\gamma}_U = (\\hat{\\gamma}_1 + \\hat{\\gamma}_2 - \\hat{\\gamma}_3) + \\sqrt{(\\tilde{\\gamma}_1 - \\hat{\\gamma}_1)^2 + (\\tilde{\\gamma}_2 - \\hat{\\gamma}_2)^2 + (\\tilde{\\gamma}_3 - \\hat{\\gamma}_3)^2}$$"),
        p(strong("Final Result = "), strong(sprintf("%.5f", res$gamma_U_hat)))
      )
    )
  })
}

# Run the Shiny application
shinyApp(ui = ui, server = server)