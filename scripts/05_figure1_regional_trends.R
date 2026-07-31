# ============================================================
# 05_figure1_regional_trends.R
#
# Purpose:
# Generate publication-quality two-panel regional resistance
# trend plots for meropenem and imipenem, with Wilson 95% CIs.
#
# Outputs:
# output/figures/Figure1_Regional_Trends.pdf
# output/figures/Figure1_Regional_Trends.png
# output/figures/Figure1_Regional_Trends.tiff
# ============================================================


# ---- Packages ------------------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "tidyr",
  "ggplot2"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)


# ---- File paths ----------------------------------------------

input_file <-
  "data/processed/carbapenem_pooled_summary.csv"

pdf_output <-
  "output/figures/Figure1_Regional_Trends.pdf"

png_output <-
  "output/figures/Figure1_Regional_Trends.png"

tiff_output <-
  "output/figures/Figure1_Regional_Trends.tiff"


# ---- Check input ---------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file,
    "\nRun scripts/02_pooled_resistance_wilson_ci.R first."
  )
}

dir.create(
  "output/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


# ---- Load data -----------------------------------------------

pooled <- read_csv(
  input_file,
  show_col_types = FALSE
)


# ---- Validate columns ----------------------------------------

required_columns <- c(
  "WHORegionName",
  "Year",
  "AntibioticName",
  "Resistance",
  "LowerCI",
  "UpperCI"
)

missing_columns <- setdiff(
  required_columns,
  names(pooled)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ---- Ordering -------------------------------------------------

region_order <- c(
  "African Region",
  "Eastern Mediterranean Region",
  "European Region",
  "Region of the Americas",
  "South-East Asia Region",
  "Western Pacific Region"
)

antibiotic_order <- c(
  "Meropenem",
  "Imipenem"
)


# ---- Colourblind-safe Okabe-Ito palette ----------------------

region_colours <- c(
  "African Region" = "#D55E00",
  "Eastern Mediterranean Region" = "#E69F00",
  "European Region" = "#0072B2",
  "Region of the Americas" = "#009E73",
  "South-East Asia Region" = "#CC79A7",
  "Western Pacific Region" = "#56B4E9"
)


# ---- Prepare plotting data -----------------------------------

plot_data <- pooled %>%
  transmute(
    WHORegionName = as.character(WHORegionName),
    Year = as.integer(Year),
    AntibioticName = as.character(AntibioticName),
    Resistance = as.numeric(Resistance),
    LowerCI = as.numeric(LowerCI),
    UpperCI = as.numeric(UpperCI)
  ) %>%
  filter(
    WHORegionName %in% region_order,
    AntibioticName %in% antibiotic_order,
    Year %in% 2020:2023
  ) %>%
  mutate(
    WHORegionName = factor(
      WHORegionName,
      levels = region_order
    ),
    AntibioticName = factor(
      AntibioticName,
      levels = antibiotic_order
    )
  ) %>%
  complete(
    WHORegionName,
    AntibioticName,
    Year = 2020:2023
  ) %>%
  arrange(
    AntibioticName,
    WHORegionName,
    Year
  )


# ---- Validation ----------------------------------------------

if (
  any(
    plot_data$Resistance < 0 |
    plot_data$Resistance > 100,
    na.rm = TRUE
  )
) {
  stop("Resistance values must be between 0 and 100.")
}

if (
  any(
    plot_data$LowerCI > plot_data$UpperCI,
    na.rm = TRUE
  )
) {
  stop("At least one lower CI exceeds the upper CI.")
}


# ---- Create Figure 1 -----------------------------------------

figure1 <- ggplot(
  plot_data,
  aes(
    x = Year,
    y = Resistance,
    colour = WHORegionName,
    group = WHORegionName
  )
) +
  geom_errorbar(
    aes(
      ymin = LowerCI,
      ymax = UpperCI
    ),
    width = 0.06,
    linewidth = 0.55,
    alpha = 0.55,
    na.rm = TRUE
  ) +
  geom_line(
    linewidth = 1.15,
    na.rm = FALSE
  ) +
  geom_point(
    size = 3.2,
    stroke = 0.3,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ AntibioticName,
    nrow = 1
  ) +
  scale_colour_manual(
    values = region_colours,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = 2020:2023,
    limits = c(2019.82, 2023.18)
  ) +
  scale_y_continuous(
    limits = c(20, 90),
    breaks = seq(20, 90, by = 10),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    x = "Year",
    y = "Resistance (%)",
    colour = "WHO region"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 13,
      face = "bold",
      margin = margin(b = 8)
    ),
    
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.6,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    panel.spacing = grid::unit(
      1.3,
      "lines"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.key.width = grid::unit(
      1.6,
      "cm"
    ),
    
    plot.margin = margin(
      t = 12,
      r = 12,
      b = 12,
      l = 12
    )
  ) +
  guides(
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.4,
        size = 3.4
      )
    )
  )


# ---- Export ---------------------------------------------------

ggsave(
  filename = pdf_output,
  plot = figure1,
  width = 12,
  height = 7.2,
  units = "in",
  device = "pdf"
)

ggsave(
  filename = png_output,
  plot = figure1,
  width = 12,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = tiff_output,
  plot = figure1,
  width = 12,
  height = 7.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ---- Confirm outputs -----------------------------------------

stopifnot(
  file.exists(pdf_output),
  file.exists(png_output),
  file.exists(tiff_output)
)

message("\nFigure 1 generated successfully.")
message("PDF: ", pdf_output)
message("PNG: ", png_output)
message("TIFF: ", tiff_output)

print(figure1)

