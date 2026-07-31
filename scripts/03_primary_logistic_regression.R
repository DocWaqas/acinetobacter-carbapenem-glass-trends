# ============================================================
# 03_primary_logistic_regression.R
#
# Purpose:
# Evaluate annual temporal trends in meropenem and imipenem
# resistance for each WHO region using binomial logistic
# regression.
#
# Model:
# cbind(Resistant, InterpretableAST - Resistant) ~ YearIndex
#
# YearIndex:
# 2020 = 0, 2021 = 1, 2022 = 2, 2023 = 3
#
# The exponentiated YearIndex coefficient is the odds ratio
# representing the annual change in the odds of resistance.
#
# Outputs:
# 1. data/processed/logistic_regression_results.csv
# 2. output/tables/Table2_Primary_Logistic_Regression.csv
# 3. output/tables/Table2_Primary_Logistic_Regression.xlsx
# ============================================================


# ---- Required packages --------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "writexl"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(dplyr)
library(readr)
library(writexl)


# ---- File paths ---------------------------------------------

input_file <-
  "data/processed/glass_carbapenem_validated.csv"

processed_output <-
  "data/processed/logistic_regression_results.csv"

table_csv_output <-
  "output/tables/Table2_Primary_Logistic_Regression.csv"

table_excel_output <-
  "output/tables/Table2_Primary_Logistic_Regression.xlsx"


# ---- Check input and output folders -------------------------

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file,
    "\nRun scripts/01_load_and_validate_data.R first."
  )
}

dir.create(
  "output/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


# ---- Load validated dataset --------------------------------

glass_carb <- read_csv(
  input_file,
  show_col_types = FALSE
)


# ---- Validate required columns ------------------------------

required_columns <- c(
  "WHORegionName",
  "Year",
  "CountryTerritoryArea",
  "AntibioticName",
  "InterpretableAST",
  "Resistant"
)

missing_columns <- setdiff(
  required_columns,
  names(glass_carb)
)

if (length(missing_columns) > 0) {
  stop(
    "The validated dataset is missing these columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ---- Prepare analysis data ----------------------------------

analysis_data <- glass_carb %>%
  mutate(
    WHORegionName = as.character(WHORegionName),
    AntibioticName = as.character(AntibioticName),
    Year = as.integer(Year),
    YearIndex = Year - 2020L,
    InterpretableAST = as.integer(InterpretableAST),
    Resistant = as.integer(Resistant),
    NonResistant = InterpretableAST - Resistant
  ) %>%
  filter(
    Year %in% 2020:2023,
    AntibioticName %in% c(
      "Meropenem",
      "Imipenem"
    ),
    !is.na(WHORegionName),
    !is.na(YearIndex),
    !is.na(InterpretableAST),
    !is.na(Resistant),
    InterpretableAST > 0,
    Resistant >= 0,
    NonResistant >= 0
  ) %>%
  arrange(
    WHORegionName,
    AntibioticName,
    Year,
    CountryTerritoryArea
  )


# ---- Validation checks --------------------------------------

if (nrow(analysis_data) == 0) {
  stop(
    "No eligible observations remained for regression."
  )
}

if (
  any(
    analysis_data$Resistant >
    analysis_data$InterpretableAST
  )
) {
  stop(
    "At least one resistant count exceeds ",
    "the interpretable AST count."
  )
}

model_groups <- analysis_data %>%
  distinct(
    WHORegionName,
    AntibioticName
  ) %>%
  arrange(
    WHORegionName,
    AntibioticName
  )

expected_groups <- 12L

if (nrow(model_groups) != expected_groups) {
  warning(
    "Expected 12 region-antibiotic groups, but found ",
    nrow(model_groups),
    ". Check whether any data are missing."
  )
}


# ---- Function to fit one regional antibiotic model ----------

fit_trend_model <- function(
    data,
    region_name,
    antibiotic_name
) {
  
  dat <- data %>%
    filter(
      WHORegionName == region_name,
      AntibioticName == antibiotic_name
    )
  
  number_rows <- nrow(dat)
  number_years <- n_distinct(dat$Year)
  number_countries <-
    n_distinct(dat$CountryTerritoryArea)
  
  if (number_rows < 2) {
    warning(
      "Skipping ",
      region_name,
      " - ",
      antibiotic_name,
      ": fewer than two observations."
    )
    
    return(NULL)
  }
  
  if (number_years < 2) {
    warning(
      "Skipping ",
      region_name,
      " - ",
      antibiotic_name,
      ": fewer than two unique years."
    )
    
    return(NULL)
  }
  
  total_tested <- sum(
    dat$InterpretableAST,
    na.rm = TRUE
  )
  
  total_resistant <- sum(
    dat$Resistant,
    na.rm = TRUE
  )
  
  model <- tryCatch(
    glm(
      cbind(
        Resistant,
        NonResistant
      ) ~ YearIndex,
      family = binomial(
        link = "logit"
      ),
      data = dat
    ),
    error = function(e) {
      warning(
        "Model failed for ",
        region_name,
        " - ",
        antibiotic_name,
        ": ",
        conditionMessage(e)
      )
      
      return(NULL)
    }
  )
  
  if (is.null(model)) {
    return(NULL)
  }
  
  coefficient_table <- summary(model)$coefficients
  
  if (
    !"YearIndex" %in%
    rownames(coefficient_table)
  ) {
    warning(
      "YearIndex coefficient not found for ",
      region_name,
      " - ",
      antibiotic_name,
      "."
    )
    
    return(NULL)
  }
  
  beta <- coefficient_table[
    "YearIndex",
    "Estimate"
  ]
  
  standard_error <- coefficient_table[
    "YearIndex",
    "Std. Error"
  ]
  
  z_value <- coefficient_table[
    "YearIndex",
    "z value"
  ]
  
  p_value <- coefficient_table[
    "YearIndex",
    "Pr(>|z|)"
  ]
  
  # Wald 95% confidence interval.
  # This matches the coefficient-based confidence intervals
  # used for the manuscript trend table.
  lower_beta <-
    beta - qnorm(0.975) * standard_error
  
  upper_beta <-
    beta + qnorm(0.975) * standard_error
  
  odds_ratio <- exp(beta)
  lower_ci <- exp(lower_beta)
  upper_ci <- exp(upper_beta)
  
  data.frame(
    Region = region_name,
    Drug = antibiotic_name,
    Rows = number_rows,
    YearsObserved = number_years,
    Countries = number_countries,
    TotalTested = total_tested,
    TotalResistant = total_resistant,
    Beta = beta,
    StandardError = standard_error,
    Z = z_value,
    OR = odds_ratio,
    LowerCI = lower_ci,
    UpperCI = upper_ci,
    P = p_value,
    NullDeviance = model$null.deviance,
    ResidualDeviance = model$deviance,
    DegreesFreedomResidual =
      model$df.residual,
    AIC = AIC(model),
    stringsAsFactors = FALSE
  )
}


# ---- Fit all region-antibiotic models ------------------------

model_results_list <- vector(
  mode = "list",
  length = nrow(model_groups)
)

for (i in seq_len(nrow(model_groups))) {
  
  region_name <-
    model_groups$WHORegionName[i]
  
  antibiotic_name <-
    model_groups$AntibioticName[i]
  
  cat(
    "Running: ",
    region_name,
    " - ",
    antibiotic_name,
    "\n",
    sep = ""
  )
  
  model_results_list[[i]] <-
    fit_trend_model(
      data = analysis_data,
      region_name = region_name,
      antibiotic_name = antibiotic_name
    )
}

model_results_list <-
  model_results_list[
    !vapply(
      model_results_list,
      is.null,
      logical(1)
    )
  ]

if (length(model_results_list) == 0) {
  stop(
    "No logistic regression models were successfully fitted."
  )
}

logistic_results <- bind_rows(
  model_results_list
)


# ---- Apply Benjamini-Hochberg adjustment --------------------

logistic_results <- logistic_results %>%
  mutate(
    AdjustedP = p.adjust(
      P,
      method = "BH"
    )
  )


# ---- Add annual percentage change in odds -------------------

logistic_results <- logistic_results %>%
  mutate(
    AnnualOddsChangePercent =
      100 * (OR - 1)
  )


# ---- Add interpretation column ------------------------------

logistic_results <- logistic_results %>%
  mutate(
    Interpretation = case_when(
      
      AdjustedP >= 0.05 ~
        "No significant temporal trend",
      
      AdjustedP < 0.05 &
        OR < 1 ~
        paste0(
          "Significant annual decrease (",
          sprintf(
            "%.1f",
            100 * (1 - OR)
          ),
          "% lower odds of resistance per year)"
        ),
      
      AdjustedP < 0.05 &
        OR > 1 ~
        paste0(
          "Significant annual increase (",
          sprintf(
            "%.1f",
            100 * (OR - 1)
          ),
          "% higher odds of resistance per year)"
        ),
      
      TRUE ~
        "No significant temporal trend"
    )
  )


# ---- Arrange rows in manuscript order -----------------------

region_order <- c(
  "African Region",
  "Eastern Mediterranean Region",
  "European Region",
  "Region of the Americas",
  "South-East Asia Region",
  "Western Pacific Region"
)

drug_order <- c(
  "Meropenem",
  "Imipenem"
)

logistic_results <- logistic_results %>%
  mutate(
    Region = factor(
      Region,
      levels = region_order
    ),
    Drug = factor(
      Drug,
      levels = drug_order
    )
  ) %>%
  arrange(
    Region,
    Drug
  ) %>%
  mutate(
    Region = as.character(Region),
    Drug = as.character(Drug)
  )


# ---- Save detailed regression results -----------------------

write_csv(
  logistic_results,
  processed_output
)


# ---- Create publication-formatted Table 2 -------------------

format_p_value <- function(x) {
  
  if (is.na(x)) {
    return(NA_character_)
  }
  
  if (x < 0.001) {
    return("<0.001")
  }
  
  sprintf("%.3f", x)
}


table2 <- logistic_results %>%
  mutate(
    `Odds Ratio (OR) per Year` =
      sprintf("%.3f", OR),
    
    `95% Confidence Interval` =
      paste0(
        sprintf("%.3f", LowerCI),
        "–",
        sprintf("%.3f", UpperCI)
      ),
    
    `p-value` =
      vapply(
        P,
        format_p_value,
        character(1)
      ),
    
    `BH-adjusted p-value` =
      vapply(
        AdjustedP,
        format_p_value,
        character(1)
      )
  ) %>%
  select(
    `WHO Region` = Region,
    Antibiotic = Drug,
    `Odds Ratio (OR) per Year`,
    `95% Confidence Interval`,
    `p-value`,
    `BH-adjusted p-value`,
    Interpretation
  )


# ---- Export Table 2 -----------------------------------------

write_csv(
  table2,
  table_csv_output
)

writexl::write_xlsx(
  table2,
  table_excel_output
)


# ---- Display completion summary -----------------------------

cat(
  "\nPrimary logistic regression completed successfully.\n"
)

cat(
  "Models fitted: ",
  nrow(logistic_results),
  "\n",
  sep = ""
)

cat(
  "Detailed results written to: ",
  processed_output,
  "\n",
  sep = ""
)

cat(
  "Table 2 CSV written to: ",
  table_csv_output,
  "\n",
  sep = ""
)

cat(
  "Table 2 Excel file written to: ",
  table_excel_output,
  "\n\n",
  sep = ""
)

print(table2)

if (nrow(logistic_results) != 12) {
  warning(
    "Expected 12 region-antibiotic models, but fitted ",
    nrow(logistic_results),
    "."
  )
}
stopifnot(nrow(table2) == 12)

rm(list = ls())
