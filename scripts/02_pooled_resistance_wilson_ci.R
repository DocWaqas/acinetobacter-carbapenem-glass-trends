# ============================================================
# 02_pooled_resistance_wilson_ci.R
#
# Purpose:
# Calculate isolate-weighted pooled carbapenem resistance
# estimates and Wilson 95% confidence intervals by WHO region,
# year, and antibiotic.
#
# Outputs:
# 1. data/processed/carbapenem_pooled_summary.csv
# 2. output/tables/Table1_Regional_Resistance_WilsonCI.csv
# 3. output/tables/Table1_Regional_Resistance_WilsonCI.xlsx
# ============================================================

# ---- Required packages --------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "binom",
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
library(binom)
library(writexl)

# ---- File paths ---------------------------------------------

input_file <- "data/processed/glass_carbapenem_validated.csv"

processed_output <- "data/processed/carbapenem_pooled_summary.csv"

table_csv_output <-
  "output/tables/Table1_Regional_Resistance_WilsonCI.csv"

table_excel_output <-
  "output/tables/Table1_Regional_Resistance_WilsonCI.xlsx"

# ---- Check input and output folders -------------------------

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
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

# ---- Create pooled regional summary -------------------------

summary_table <- glass_carb %>%
  group_by(
    WHORegionName,
    Year,
    AntibioticName
  ) %>%
  summarise(
    TotalTested = sum(
      InterpretableAST,
      na.rm = TRUE
    ),
    TotalResistant = sum(
      Resistant,
      na.rm = TRUE
    ),
    Countries = n_distinct(
      CountryTerritoryArea
    ),
    .groups = "drop"
  ) %>%
  arrange(
    WHORegionName,
    Year,
    AntibioticName
  )

# ---- Validate aggregated counts -----------------------------

if (any(summary_table$TotalTested <= 0)) {
  stop(
    "At least one pooled group has TotalTested <= 0."
  )
}

if (
  any(
    summary_table$TotalResistant >
    summary_table$TotalTested
  )
) {
  stop(
    "At least one pooled group has resistant counts ",
    "greater than tested counts."
  )
}

# ---- Wilson confidence intervals ----------------------------

wilson_ci <- binom::binom.confint(
  x = summary_table$TotalResistant,
  n = summary_table$TotalTested,
  methods = "wilson"
)

if (nrow(wilson_ci) != nrow(summary_table)) {
  stop(
    "Wilson CI output does not match the number of ",
    "rows in the pooled summary."
  )
}

# ---- Add resistance estimates and confidence limits ---------

summary_table <- summary_table %>%
  mutate(
    Resistance = 100 *
      TotalResistant /
      TotalTested,
    LowerCI = 100 *
      wilson_ci$lower,
    UpperCI = 100 *
      wilson_ci$upper
  )

# ---- Save detailed processed summary ------------------------

write_csv(
  summary_table,
  processed_output
)

# ---- Create publication-formatted Table 1 -------------------

table1 <- summary_table %>%
  mutate(
    `Resistance % (95% CI)` = paste0(
      sprintf("%.1f", Resistance),
      " (",
      sprintf("%.1f", LowerCI),
      "–",
      sprintf("%.1f", UpperCI),
      ")"
    )
  ) %>%
  select(
    `WHO Region` = WHORegionName,
    Year,
    Antibiotic = AntibioticName,
    Countries,
    `Total Tested` = TotalTested,
    `Total Resistant` = TotalResistant,
    `Resistance % (95% CI)`
  ) %>%
  arrange(
    `WHO Region`,
    Year,
    Antibiotic
  )

# ---- Export Table 1 -----------------------------------------

write_csv(
  table1,
  table_csv_output
)

writexl::write_xlsx(
  table1,
  table_excel_output
)

# ---- Display results ----------------------------------------

cat("\nPooled resistance analysis completed successfully.\n")
cat(
  "Number of pooled region-year-antibiotic rows: ",
  nrow(summary_table),
  "\n",
  sep = ""
)
cat(
  "Processed summary written to: ",
  processed_output,
  "\n",
  sep = ""
)
cat(
  "Table 1 CSV written to: ",
  table_csv_output,
  "\n",
  sep = ""
)
cat(
  "Table 1 Excel file written to: ",
  table_excel_output,
  "\n\n",
  sep = ""
)

print(
  table1,
  n = Inf,
  width = Inf
)
