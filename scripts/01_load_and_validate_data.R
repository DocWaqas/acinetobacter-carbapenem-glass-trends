# ============================================================
# 01_load_and_validate_data.R
#
# Purpose:
# Load and validate the processed WHO GLASS dataset used for
# the carbapenem-resistance analysis.
# ============================================================

# ---- Packages -----------------------------------------------

required_packages <- c("dplyr", "readr")

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(dplyr)
library(readr)

# ---- Paths --------------------------------------------------

input_file <- "data/processed/glass_combined.csv"
output_file <- "data/processed/glass_carbapenem_validated.csv"

if (!file.exists(input_file)) {
  stop(
    "File not found: ", input_file,
    "\nOpen the RStudio project from the repository root."
  )
}

# ---- Load data ----------------------------------------------

glass <- read_csv(
  input_file,
  show_col_types = FALSE
)

# ---- Required columns ---------------------------------------

required_columns <- c(
  "WHORegionName",
  "Year",
  "Iso3",
  "CountryTerritoryArea",
  "Specimen",
  "PathogenName",
  "AntibioticName",
  "InterpretableAST",
  "Resistant",
  "ResistancePercentage",
  "SourceFile"
)

missing_columns <- setdiff(required_columns, names(glass))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ---- Standardize variable types -----------------------------

glass <- glass %>%
  mutate(
    WHORegionName = as.character(WHORegionName),
    Year = as.integer(Year),
    Iso3 = as.character(Iso3),
    CountryTerritoryArea = as.character(CountryTerritoryArea),
    Specimen = as.character(Specimen),
    PathogenName = as.character(PathogenName),
    AntibioticName = as.character(AntibioticName),
    InterpretableAST = as.integer(InterpretableAST),
    Resistant = as.integer(Resistant),
    ResistancePercentage = as.numeric(ResistancePercentage),
    SourceFile = as.character(SourceFile)
  )

# ---- Restrict to manuscript scope ---------------------------

glass_carb <- glass %>%
  filter(
    Year %in% 2020:2023,
    Specimen == "Bloodstream",
    PathogenName == "Acinetobacter spp.",
    AntibioticName %in% c("Meropenem", "Imipenem"),
    !is.na(WHORegionName),
    !is.na(CountryTerritoryArea),
    !is.na(InterpretableAST),
    !is.na(Resistant),
    InterpretableAST > 0,
    Resistant >= 0,
    Resistant <= InterpretableAST
  ) %>%
  arrange(
    WHORegionName,
    Year,
    AntibioticName,
    CountryTerritoryArea
  )

# ---- Validation checks --------------------------------------

if (nrow(glass_carb) == 0) {
  stop("No eligible observations remained after filtering.")
}

if (any(glass_carb$Resistant > glass_carb$InterpretableAST)) {
  stop("At least one resistant count exceeds the tested count.")
}

# ---- Save validated dataset ---------------------------------

write_csv(
  glass_carb,
  output_file
)

# ---- Summary ------------------------------------------------

cat("\nValidation completed successfully.\n")
cat("Rows in combined dataset: ", nrow(glass), "\n", sep = "")
cat("Rows in validated carbapenem dataset: ",
    nrow(glass_carb), "\n", sep = "")
cat("Saved to: ", output_file, "\n\n", sep = "")

validation_summary <- glass_carb %>%
  summarise(
    Observations = n(),
    Regions = n_distinct(WHORegionName),
    Countries = n_distinct(CountryTerritoryArea),
    Years = n_distinct(Year),
    Antibiotics = n_distinct(AntibioticName)
  )

print(validation_summary)

reporting_summary <- glass_carb %>%
  count(
    WHORegionName,
    Year,
    AntibioticName,
    name = "Rows"
  ) %>%
  arrange(
    WHORegionName,
    Year,
    AntibioticName
  )

print(reporting_summary, n = Inf)

rm(list = ls())

source("scripts/01_load_and_validate_data.R")
file.exists("data/processed/glass_combined.csv")

list.files("data/processed")

getwd()

list.files("data/processed")

source("scripts/01_load_and_validate_data.R")
