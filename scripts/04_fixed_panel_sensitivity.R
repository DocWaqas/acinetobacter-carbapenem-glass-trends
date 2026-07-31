# ============================================================
# 04_fixed_panel_sensitivity.R
#
# Purpose:
# Repeat the regional temporal trend analysis using only
# countries that reported continuously in every year from
# 2020 to 2023 for each antibiotic.
#
# Outputs:
# 1. data/processed/fixed_panel_regression_results.csv
# 2. data/processed/fixed_panel_countries.csv
# 3. output/tables/Supplementary_Table_S2_Fixed_Panel.csv
# 4. output/tables/Supplementary_Table_S2_Fixed_Panel.xlsx
# ============================================================


# ---- Packages ------------------------------------------------

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


# ---- File paths ----------------------------------------------

input_file <-
  "data/processed/glass_carbapenem_validated.csv"

results_output <-
  "data/processed/fixed_panel_regression_results.csv"

countries_output <-
  "data/processed/fixed_panel_countries.csv"

table_csv_output <-
  "output/tables/Supplementary_Table_S2_Fixed_Panel.csv"

table_excel_output <-
  "output/tables/Supplementary_Table_S2_Fixed_Panel.xlsx"


# ---- Check files and folders ---------------------------------

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file,
    "\nRun scripts/01_load_and_validate_data.R first."
  )
}

dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "output/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


# ---- Load validated data -------------------------------------

glass_carb <- read_csv(
  input_file,
  show_col_types = FALSE
)


# ---- Validate required columns -------------------------------

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
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ---- Prepare analysis data -----------------------------------

analysis_data <- glass_carb %>%
  transmute(
    WHORegionName = as.character(WHORegionName),
    Year = as.integer(Year),
    CountryTerritoryArea =
      as.character(CountryTerritoryArea),
    AntibioticName =
      as.character(AntibioticName),
    InterpretableAST =
      as.integer(InterpretableAST),
    Resistant =
      as.integer(Resistant),
    NonResistant =
      as.integer(InterpretableAST) -
      as.integer(Resistant),
    YearIndex =
      as.integer(Year) - 2020L
  ) %>%
  filter(
    Year %in% 2020:2023,
    AntibioticName %in% c(
      "Meropenem",
      "Imipenem"
    ),
    !is.na(WHORegionName),
    !is.na(CountryTerritoryArea),
    !is.na(InterpretableAST),
    !is.na(Resistant),
    InterpretableAST > 0,
    Resistant >= 0,
    NonResistant >= 0
  )


# ---- Identify continuously reporting countries --------------

country_year_counts <- analysis_data %>%
  distinct(
    WHORegionName,
    AntibioticName,
    CountryTerritoryArea,
    Year
  ) %>%
  count(
    WHORegionName,
    AntibioticName,
    CountryTerritoryArea,
    name = "YearsReported"
  )

fixed_countries <- country_year_counts %>%
  filter(
    YearsReported == 4
  ) %>%
  arrange(
    WHORegionName,
    AntibioticName,
    CountryTerritoryArea
  )

if (nrow(fixed_countries) == 0) {
  stop(
    "No countries reported continuously across all four years."
  )
}

write_csv(
  fixed_countries,
  countries_output
)


# ---- Create fixed-panel dataset ------------------------------

fixed_panel <- analysis_data %>%
  semi_join(
    fixed_countries,
    by = c(
      "WHORegionName",
      "AntibioticName",
      "CountryTerritoryArea"
    )
  ) %>%
  arrange(
    WHORegionName,
    AntibioticName,
    CountryTerritoryArea,
    Year
  )

if (nrow(fixed_panel) == 0) {
  stop(
    "The fixed-panel dataset contains no observations."
  )
}


# ---- Validate four-year continuity ---------------------------

continuity_check <- fixed_panel %>%
  distinct(
    WHORegionName,
    AntibioticName,
    CountryTerritoryArea,
    Year
  ) %>%
  count(
    WHORegionName,
    AntibioticName,
    CountryTerritoryArea,
    name = "YearsObserved"
  )

if (any(continuity_check$YearsObserved != 4)) {
  stop(
    "At least one fixed-panel country does not have four years ",
    "of observations."
  )
}


# ---- Display panel composition -------------------------------

panel_summary <- fixed_panel %>%
  summarise(
    Observations = n(),
    UniqueCountries =
      n_distinct(CountryTerritoryArea),
    RegionAntibioticPanels =
      n_distinct(
        paste(
          WHORegionName,
          AntibioticName
        )
      )
  )

message("\nFixed-panel dataset summary:")
print(panel_summary)

panel_counts <- fixed_panel %>%
  distinct(
    WHORegionName,
    AntibioticName,
    CountryTerritoryArea
  ) %>%
  count(
    WHORegionName,
    AntibioticName,
    name = "Countries"
  ) %>%
  arrange(
    WHORegionName,
    AntibioticName
  )

print(panel_counts)


# ---- Model-fitting function ----------------------------------

fit_fixed_panel_model <- function(
    data,
    region_name,
    antibiotic_name
) {
  
  dat <- data %>%
    filter(
      WHORegionName == region_name,
      AntibioticName == antibiotic_name
    )
  
  if (nrow(dat) == 0) {
    return(NULL)
  }
  
  if (n_distinct(dat$Year) != 4) {
    warning(
      "Skipping ",
      region_name,
      " - ",
      antibiotic_name,
      ": not all four years are represented."
    )
    return(NULL)
  }
  
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
  
  if (!"YearIndex" %in% rownames(coefficient_table)) {
    warning(
      "YearIndex coefficient missing for ",
      region_name,
      " - ",
      antibiotic_name
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
  
  p_value <- coefficient_table[
    "YearIndex",
    "Pr(>|z|)"
  ]
  
  z_critical <- qnorm(0.975)
  
  lower_beta <-
    beta - z_critical * standard_error
  
  upper_beta <-
    beta + z_critical * standard_error
  
  tibble(
    Region = region_name,
    Drug = antibiotic_name,
    Rows = nrow(dat),
    Countries =
      n_distinct(dat$CountryTerritoryArea),
    TotalTested =
      sum(dat$InterpretableAST),
    TotalResistant =
      sum(dat$Resistant),
    Beta = beta,
    StandardError = standard_error,
    OR = exp(beta),
    LowerCI = exp(lower_beta),
    UpperCI = exp(upper_beta),
    P = p_value,
    ResidualDeviance = model$deviance,
    ResidualDF = model$df.residual,
    AIC = AIC(model)
  )
}


# ---- Fit all fixed-panel models ------------------------------

model_groups <- fixed_panel %>%
  distinct(
    WHORegionName,
    AntibioticName
  ) %>%
  arrange(
    WHORegionName,
    AntibioticName
  )

results_list <- vector(
  mode = "list",
  length = nrow(model_groups)
)

for (i in seq_len(nrow(model_groups))) {
  
  region_name <-
    model_groups$WHORegionName[i]
  
  antibiotic_name <-
    model_groups$AntibioticName[i]
  
  message(
    "Running fixed-panel model: ",
    region_name,
    " - ",
    antibiotic_name
  )
  
  results_list[[i]] <-
    fit_fixed_panel_model(
      data = fixed_panel,
      region_name = region_name,
      antibiotic_name = antibiotic_name
    )
}

results_list <- results_list[
  !vapply(
    results_list,
    is.null,
    logical(1)
  )
]

if (length(results_list) == 0) {
  stop(
    "No fixed-panel models were successfully fitted."
  )
}

logistic_results_fixed <- bind_rows(
  results_list
)


# ---- Apply BH adjustment -------------------------------------

logistic_results_fixed <- logistic_results_fixed %>%
  mutate(
    AdjustedP = p.adjust(
      P,
      method = "BH"
    )
  )


# ---- Add interpretations -------------------------------------

logistic_results_fixed <- logistic_results_fixed %>%
  mutate(
    Interpretation = case_when(
      AdjustedP >= 0.05 ~
        "No significant temporal trend",
      
      AdjustedP < 0.05 & OR < 1 ~
        paste0(
          "Significant annual decrease (",
          sprintf(
            "%.1f",
            100 * (1 - OR)
          ),
          "% lower odds of resistance per year)"
        ),
      
      AdjustedP < 0.05 & OR > 1 ~
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


# ---- Arrange manuscript order --------------------------------

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

logistic_results_fixed <- logistic_results_fixed %>%
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


# ---- Validate expected model count ---------------------------

if (nrow(logistic_results_fixed) != 10) {
  warning(
    "Expected 10 fixed-panel region-antibiotic models, but fitted ",
    nrow(logistic_results_fixed),
    "."
  )
}


# ---- Export detailed results ---------------------------------

write_csv(
  logistic_results_fixed,
  results_output
)


# ---- Format publication table --------------------------------

format_p_value <- function(x) {
  
  if (is.na(x)) {
    return(NA_character_)
  }
  
  if (x < 0.001) {
    return("<0.001")
  }
  
  sprintf("%.3f", x)
}

supp_table_s2 <- logistic_results_fixed %>%
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
    Countries,
    `Odds Ratio (OR) per Year`,
    `95% Confidence Interval`,
    `p-value`,
    `BH-adjusted p-value`,
    Interpretation
  )


# ---- Export Supplementary Table S2 ---------------------------

write_csv(
  supp_table_s2,
  table_csv_output
)

write_xlsx(
  supp_table_s2,
  table_excel_output
)


# ---- Final validation ----------------------------------------

stopifnot(
  nrow(logistic_results_fixed) ==
    nrow(supp_table_s2)
)

message(
  "\nFixed-panel sensitivity analysis completed successfully."
)

message(
  "Fixed-panel observations: ",
  nrow(fixed_panel)
)

message(
  "Models fitted: ",
  nrow(logistic_results_fixed)
)

message(
  "Fixed-panel country list written to: ",
  countries_output
)

message(
  "Detailed results written to: ",
  results_output
)

message(
  "Supplementary Table S2 CSV written to: ",
  table_csv_output
)

message(
  "Supplementary Table S2 Excel written to: ",
  table_excel_output
)

print(supp_table_s2)